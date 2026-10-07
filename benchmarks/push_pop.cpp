#include "util/benchmark.hpp"
#include "util/memory_stats.hpp"
#include "util/thread_coordination.hpp"
#include "wrapper/selector.hpp"

#include <cxxopts.hpp>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <random>
#include <utility>
#include <vector>

using key_type = unsigned long;
using value_type = unsigned long;

using clock_type = benchmark::clock_type;

using pq_type = PQ<true, key_type, value_type>;
using handle_type = pq_type::handle_type;

struct Settings {
    benchmark::BaseSettings base_settings{};
    pq_type::settings_type pq_settings;
    long long prefill_per_thread = 1 << 20;
    long long elements_per_thread = 1 << 14;
    long long batch_size = 1 << 12;
    int seed = 1;
    benchmark::Papi papi;

    void register_cmd_options(cxxopts::Options& cmd) {
        base_settings.register_cmd_options(cmd);
        pq_settings.register_cmd_options(cmd);
        papi.register_cmd_options(cmd);
        cmd.add_options()
            // clang-format off
            ("p,prefill", "Prefill per thread", cxxopts::value<long long>(prefill_per_thread), "NUMBER")
            ("n,elements", "Number of elements per thread", cxxopts::value<long long>(elements_per_thread), "NUMBER")
            ("batch-size", "Batch size", cxxopts::value<long long>(batch_size), "NUMBER")
            ("s,seed", "Initial seed", cxxopts::value<int>(seed), "NUMBER")
            // clang-format on
            ;
    }

    bool validate() const {
        if (!base_settings.validate()) {
            return false;
        }
        if (!pq_settings.validate()) {
            return false;
        }
        if (prefill_per_thread < 0) {
            std::cerr << "Error: Prefill must be nonnegative\n";
            return false;
        }
        if (elements_per_thread < 0) {
            std::cerr << "Error: Number of elements must be nonnegative\n";
            return false;
        }
        if (batch_size <= 0) {
            std::cerr << "Error: batch size must be greater than 0\n";
            return false;
        }
        if (seed <= 0) {
            std::cerr << "Error: Seed must be greater than 0\n";
            return false;
        }
        if (!papi.validate()) {
            return false;
        }
        return true;
    }

    void write_human_readable(std::ostream& out) const {
        base_settings.write_human_readable(out);
        pq_settings.write_human_readable(out);
        out << "Prefill per thread: " << prefill_per_thread << '\n';
        out << "Elements per thread: " << elements_per_thread << '\n';
        out << "Batch size: " << batch_size << '\n';
        out << "Seed: " << seed << '\n';
        papi.write_human_readable(out);
    }

    void write_json(json::Object& obj) const {
        base_settings.write_json(obj);
        obj.object("pq", [this](json::Object& o) { pq_settings.write_json(o); });
        obj.entry("prefill_per_thread", prefill_per_thread);
        obj.entry("elements_per_thread", elements_per_thread);
        obj.entry("batch_size", batch_size);
        obj.entry("seed", seed);
        papi.write_json(obj);
    }
};

struct ThreadData {
    benchmark::Interval push_interval{};
    benchmark::Interval pop_interval{};
    long long push_count{0};
    long long pop_count{0};
    long long failed_pop_count{0};
    int event_set = -1;
    std::vector<long long> push_papi_event_counter{};
    std::vector<long long> pop_papi_event_counter{};

    void write_json(json::Object& obj) const {
        obj.entry("pushes", push_count);
        obj.entry("pops", pop_count);
        obj.entry("failed_pops", failed_pop_count);
        obj.array("push_papi_event_counter", push_papi_event_counter);
        obj.array("pop_papi_event_counter", pop_papi_event_counter);
    }
};

struct SharedData {
    std::vector<key_type> keys;
    std::atomic_llong counter{0};
    memory_stats::Snapshot memory_start;
    memory_stats::Snapshot memory_after_push;
    memory_stats::Snapshot memory_end;
    std::vector<ThreadData> thread_data;
};

std::vector<benchmark::Interval> intervals(SharedData const& data, benchmark::Interval ThreadData::* phase) {
    std::vector<benchmark::Interval> result;
    for (auto const& t : data.thread_data) {
        result.push_back(t.*phase);
    }
    return result;
}

void write_result_json(Settings const& settings, SharedData const& data, std::ostream& out) {
    {
        json::Object root{out};
        root.object("settings", [&settings](json::Object& obj) { settings.write_json(obj); });
        root.object("results", [&data](json::Object& results) {
            benchmark::write_timing(results, "push_", intervals(data, &ThreadData::push_interval));
            benchmark::write_timing(results, "pop_", intervals(data, &ThreadData::pop_interval));
            results.object("memory", [&data](json::Object& memory) {
                memory_stats::write_json(
                    memory,
                    {{"start", data.memory_start}, {"after_push", data.memory_after_push}, {"end", data.memory_end}});
            });
            results.array("thread_data", data.thread_data.begin(), data.thread_data.end(),
                          [](std::ostream& out, ThreadData const& thread_data) {
                              json::Object obj{out};
                              thread_data.write_json(obj);
                          });
        });
    }
    out << '\n';
}

class Context : public thread_coordination::Context {
    handle_type handle_;
    ThreadData thread_data_;
    SharedData* shared_data_;
    Settings const* settings_;

   public:
    explicit Context(thread_coordination::Context ctx, handle_type handle, SharedData& shared_data,
                     Settings const& settings)
        : thread_coordination::Context{std::move(ctx)},
          handle_{std::move(handle)},
          shared_data_{&shared_data},
          settings_{&settings} {
    }

    void push(std::pair<key_type, value_type> const& e) {
        handle_.push(e);
    }

    auto try_pop() {
        return handle_.try_pop();
    }

    ThreadData& thread_data() noexcept {
        return thread_data_;
    }
    [[nodiscard]] ThreadData const& thread_data() const noexcept {
        return thread_data_;
    }
    SharedData& shared_data() noexcept {
        return *shared_data_;
    }
    [[nodiscard]] SharedData const& shared_data() const noexcept {
        return *shared_data_;
    }

    [[nodiscard]] Settings const& settings() const noexcept {
        return *settings_;
    }
};

[[gnu::noinline]] void push(Context& context) {
    auto offset =
        static_cast<value_type>(context.settings().base_settings.num_threads * context.settings().prefill_per_thread);
    long long max = context.settings().elements_per_thread * context.settings().base_settings.num_threads;
    if (context.id() == 0) {
        context.shared_data().memory_start = memory_stats::Snapshot::take();
    }
    context.synchronize();
    context.settings().papi.start(context.thread_data().event_set);
    context.spin_synchronize();
    context.thread_data().push_interval.start = clock_type::now();
    while (true) {
        auto start = context.shared_data().counter.fetch_add(context.settings().batch_size, std::memory_order_relaxed);
        if (start >= max) {
            break;
        }
        auto end = std::min(start + context.settings().batch_size, max);
        for (auto i = start; i < end; ++i) {
            context.push(
                {context.shared_data().keys[static_cast<std::size_t>(i)], offset + static_cast<value_type>(i)});
        }
        context.thread_data().push_count += end - start;
    }
    context.thread_data().push_interval.end = clock_type::now();
    context.settings().papi.stop(context.thread_data().event_set, context.thread_data().push_papi_event_counter);
    context.synchronize();
    if (context.id() == 0) {
        context.shared_data().memory_after_push = memory_stats::Snapshot::take();
    }
}

[[gnu::noinline]] void pop(Context& context) {
    if (context.id() == 0) {
        context.shared_data().counter.store(0, std::memory_order_relaxed);
    }
    auto max = context.settings().elements_per_thread * context.settings().base_settings.num_threads;
    context.synchronize();
    context.settings().papi.start(context.thread_data().event_set);
    context.spin_synchronize();
    context.thread_data().pop_interval.start = clock_type::now();
    while (true) {
        long long deletions{0};
        while (context.try_pop()) {
            ++deletions;
        }
        if (deletions == 0) {
            auto current = context.shared_data().counter.load(std::memory_order_relaxed);
            if (current >= max) {
                break;
            }
        } else {
            context.thread_data().pop_count += deletions;
            auto current = context.shared_data().counter.fetch_add(deletions, std::memory_order_relaxed) + deletions;
            if (current >= max) {
                break;
            }
        }
        ++context.thread_data().failed_pop_count;
    }
    context.thread_data().pop_interval.end = clock_type::now();
    context.settings().papi.stop(context.thread_data().event_set, context.thread_data().pop_papi_event_counter);
    context.synchronize();
    if (context.id() == 0) {
        context.shared_data().memory_end = memory_stats::Snapshot::take();
    }
}

void benchmark_thread(Context context) {
    context.thread_data().event_set = context.settings().papi.create_event_set();
    std::vector<key_type> prefill(static_cast<std::size_t>(context.settings().prefill_per_thread));
    if (context.id() == 0) {
        std::clog << "Preparing...\n";
    }
    std::seed_seq seed{context.settings().seed, context.id()};
    std::default_random_engine rng(seed);
    context.synchronize();
    std::generate(prefill.begin(), prefill.end(),
                  [&rng, min = 1UL,
                   max = static_cast<key_type>(context.settings().elements_per_thread *
                                               context.settings().base_settings.num_threads)]() {
                      return std::uniform_int_distribution<key_type>(min, max)(rng);
                  });
    std::generate_n(context.shared_data().keys.begin() + context.id() * context.settings().elements_per_thread,
                    context.settings().elements_per_thread,
                    [&rng, min = 1UL,
                     max = static_cast<key_type>(context.settings().elements_per_thread *
                                                 context.settings().base_settings.num_threads)]() {
                        return std::uniform_int_distribution<key_type>(min, max)(rng);
                    });
    context.synchronize();
    if (context.id() == 0) {
        std::clog << "Prefilling...\n";
    }
    context.synchronize();
    for (auto i = 0LL; i < context.settings().prefill_per_thread; ++i) {
        context.push({prefill[static_cast<std::size_t>(i)],
                      static_cast<value_type>(context.id() * context.settings().prefill_per_thread + i)});
    }
    context.synchronize();
    if (context.id() == 0) {
        std::clog << "Pushing...\n";
    }
    push(context);
    if (context.id() == 0) {
        std::clog << "Popping...\n";
    }
    pop(context);
    context.shared_data().thread_data[static_cast<std::size_t>(context.id())] = std::move(context.thread_data());
}

void run_benchmark(Settings const& settings) {
    auto num_threads = settings.base_settings.num_threads;
    SharedData shared_data;
    shared_data.keys.resize(static_cast<std::size_t>(num_threads * settings.elements_per_thread));
    shared_data.thread_data.resize(static_cast<std::size_t>(num_threads));

    auto pq =
        pq_type(num_threads,
                static_cast<std::size_t>(num_threads * (settings.elements_per_thread + settings.prefill_per_thread)),
                settings.pq_settings);

    thread_coordination::dispatch(settings.base_settings.affinity, num_threads, [&](auto ctx) {
        benchmark_thread(Context(std::move(ctx), pq.get_handle(), shared_data, settings));
    });

    std::clog << "Done\n";
    std::clog << '\n';
    std::clog << "= Results =\n";
    std::clog << "Push time: " << std::fixed << std::setprecision(3)
              << benchmark::seconds(intervals(shared_data, &ThreadData::push_interval)) << " s\n";
    std::clog << "Pop time: " << std::fixed << std::setprecision(3)
              << benchmark::seconds(intervals(shared_data, &ThreadData::pop_interval)) << " s\n";
    write_result_json(settings, shared_data, std::cout);
}

int main(int argc, char* argv[]) {
    return benchmark::run<pq_type, Settings>(argc, argv, run_benchmark);
}
