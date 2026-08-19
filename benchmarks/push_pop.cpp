#include "util/benchmark.hpp"
#include "util/selector.hpp"
#include "util/thread_coordination.hpp"

#include <cxxopts.hpp>

#ifdef WITH_PAPI
#include <papi.h>
#include <pthread.h>
#endif
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

using pq_type = PQ<true, key_type, value_type>;
using handle_type = pq_type::handle_type;

struct Settings : benchmark::CommonSettings {
    long long prefill_per_thread = 1 << 20;
    long long elements_per_thread = 1 << 14;
    long long batch_size = 1 << 12;
    int seed = 1;
#ifdef WITH_PAPI
    std::vector<std::string> papi_events;
#endif
    pq_type::settings_type pq_settings;
};

void register_cmd_options(Settings& settings, cxxopts::Options& cmd) {
    settings.CommonSettings::register_cmd_options(cmd);
    cmd.add_options()
        // clang-format off
            ("p,prefill", "Prefill per thread", cxxopts::value<long long>(settings.prefill_per_thread), "NUMBER")
            ("n,elements", "Number of elements per thread", cxxopts::value<long long>(settings.elements_per_thread), "NUMBER")
            ("batch-size", "Batch size", cxxopts::value<long long>(settings.batch_size), "NUMBER")
            ("s,seed", "Initial seed", cxxopts::value<int>(settings.seed), "NUMBER")
#ifdef WITH_PAPI
            ("r,count-event", "Papi event to count", cxxopts::value<std::vector<std::string>>(settings.papi_events), "STRING")
#endif
        // clang-format on
        ;
    settings.pq_settings.register_cmd_options(cmd);
}

bool validate_settings(Settings const& settings) {
    if (!settings.CommonSettings::validate()) {
        return false;
    }
    if (settings.prefill_per_thread < 0) {
        std::cerr << "Error: Prefill must be nonnegative\n";
        return false;
    }
    if (settings.elements_per_thread < 0) {
        std::cerr << "Error: Number of elements must be nonnegative\n";
        return false;
    }
    if (settings.batch_size <= 0) {
        std::cerr << "Error: batch size must be greater than 0\n";
        return false;
    }
    if (settings.elements_per_thread % settings.batch_size != 0) {
        std::cerr << "Error: Number of elements must be divisible by batch size\n";
        return false;
    }
    if (settings.seed <= 0) {
        std::cerr << "Error: Seed must be greater than 0\n";
        return false;
    }
#ifdef WITH_PAPI
    if (!settings.papi_events.empty()) {
        if (int ret = PAPI_library_init(PAPI_VER_CURRENT); ret != PAPI_VER_CURRENT) {
            std::cerr << "Error: Failed to initialize PAPI library\n";
            return false;
        }
        if (int ret = PAPI_thread_init(pthread_self); ret != PAPI_OK) {
            std::cerr << "Error: Failed to initialize PAPI thread support\n";
            return false;
        }
        for (auto const& name : settings.papi_events) {
            if (PAPI_query_named_event(name.c_str()) != PAPI_OK) {
                std::cerr << "Error: PAPI event '" << name << "' not available\n";
                return false;
            }
        }
    }
#endif
    return settings.pq_settings.validate();
}

void write_settings_human_readable(Settings const& settings, std::ostream& out) {
    settings.CommonSettings::write_human_readable(out);
    out << "Prefill per thread: " << settings.prefill_per_thread << '\n';
    out << "Elements per thread: " << settings.elements_per_thread << '\n';
    out << "Batch size: " << settings.batch_size << '\n';
    out << "Seed: " << settings.seed << '\n';
#ifdef WITH_PAPI
    out << "PAPI events:";
    if (settings.papi_events.empty()) {
        out << " None\n";
    } else {
        for (auto const& e : settings.papi_events) {
            out << ' ' << e;
        }
        out << '\n';
    }
#endif
    settings.pq_settings.write_human_readable(out);
}


void write_settings_json(Settings const& settings, json::Object& obj) {
    settings.CommonSettings::write_json(obj);
    obj.entry("prefill_per_thread", settings.prefill_per_thread);
    obj.entry("elements_per_thread", settings.elements_per_thread);
    obj.entry("batch_size", settings.batch_size);
    obj.entry("seed", settings.seed);
#ifdef WITH_PAPI
    obj.array("papi_events", settings.papi_events);
#endif
    obj.raw("pq", [&settings](std::ostream& out) { settings.pq_settings.write_json(out); });
}

struct ThreadData {
    long long push_count{0};
    long long pop_count{0};
    long long failed_pop_count{0};
#ifdef WITH_PAPI
    int event_set = PAPI_NULL;
    std::vector<long long> push_papi_event_counter{};
    std::vector<long long> pop_papi_event_counter{};
#endif
};

void write_thread_data_json(ThreadData const& data, json::Object& obj) {
    obj.entry("pushes", data.push_count);
    obj.entry("pops", data.pop_count);
    obj.entry("failed_pops", data.failed_pop_count);
#ifdef WITH_PAPI
    obj.array("push_papi_event_counter", data.push_papi_event_counter);
    obj.array("pop_papi_event_counter", data.pop_papi_event_counter);
#endif
}

struct SharedData {
    std::vector<key_type> keys;
    std::atomic_llong counter{0};
    std::chrono::high_resolution_clock::time_point start_time;
    std::chrono::nanoseconds push_time{0};
    std::chrono::nanoseconds pop_time{0};
    std::vector<ThreadData> thread_data;
};

void write_result_json(Settings const& settings, SharedData const& data, std::ostream& out) {
    {
        json::Object root{out};
        root.object("settings", [&settings](json::Object& obj) { write_settings_json(settings, obj); });
        root.object("results", [&data](json::Object& results) {
            results.entry("push_time_ns", data.push_time.count());
            results.entry("pop_time_ns", data.pop_time.count());
            results.array("thread_data", data.thread_data.begin(), data.thread_data.end(),
                          [](std::ostream& out, ThreadData const& thread_data) {
                              json::Object obj{out};
                              write_thread_data_json(thread_data, obj);
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
    auto offset = static_cast<value_type>(context.settings().num_threads * context.settings().prefill_per_thread);
    long long max = context.settings().elements_per_thread * context.settings().num_threads;
    context.synchronize();
#ifdef WITH_PAPI
    if (!context.settings().papi_events.empty()) {
        if (int ret = PAPI_start(context.thread_data().event_set); ret != PAPI_OK) {
            context.write(std::cerr) << "Failed to start performance counters\n";
        }
        context.synchronize();
    }
#endif
    if (context.id() == 0) {
        context.shared_data().start_time = std::chrono::high_resolution_clock::now();
    }
    context.synchronize();
    while (true) {
        auto start = context.shared_data().counter.fetch_add(context.settings().batch_size, std::memory_order_relaxed);
        if (start >= max) {
            break;
        }
        for (auto i = start; i < start + context.settings().batch_size; ++i) {
            context.push(
                {context.shared_data().keys[static_cast<std::size_t>(i)], offset + static_cast<value_type>(i)});
        }
        context.thread_data().push_count += context.settings().batch_size;
    }
    context.synchronize();
    if (context.id() == 0) {
        auto end_time = std::chrono::high_resolution_clock::now();
        context.shared_data().push_time = end_time - context.shared_data().start_time;
    }
#ifdef WITH_PAPI
    if (!context.settings().papi_events.empty()) {
        if (int ret = PAPI_stop(context.thread_data().event_set, context.thread_data().push_papi_event_counter.data());
            ret != PAPI_OK) {
            context.write(std::cerr) << "Failed to stop performance counters\n";
        }
    }
#endif
}

[[gnu::noinline]] void pop(Context& context) {
    context.synchronize();
#ifdef WITH_PAPI
    if (!context.settings().papi_events.empty()) {
        if (int ret = PAPI_start(context.thread_data().event_set); ret != PAPI_OK) {
            context.write(std::cerr) << "Failed to start performance counters\n";
        }
        context.synchronize();
    }
#endif
    if (context.id() == 0) {
        context.shared_data().counter.store(0, std::memory_order_relaxed);
        context.shared_data().start_time = std::chrono::high_resolution_clock::now();
    }
    auto max = context.settings().elements_per_thread * context.settings().num_threads;
    context.synchronize();
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
    context.synchronize();
    if (context.id() == 0) {
        auto end_time = std::chrono::high_resolution_clock::now();
        context.shared_data().pop_time = end_time - context.shared_data().start_time;
    }
#ifdef WITH_PAPI
    if (!context.settings().papi_events.empty()) {
        if (int ret = PAPI_stop(context.thread_data().event_set, context.thread_data().pop_papi_event_counter.data());
            ret != PAPI_OK) {
            context.write(std::cerr) << "Failed to stop performance counters\n";
        }
    }
#endif
}

#ifdef WITH_PAPI
int prepare_papi(Settings const& settings) {
    if (int ret = PAPI_register_thread(); ret != PAPI_OK) {
        throw std::runtime_error{"Failed to register thread for PAPI"};
    }
    int event_set = PAPI_NULL;
    if (int ret = PAPI_create_eventset(&event_set); ret != PAPI_OK) {
        throw std::runtime_error{"Failed to create PAPI event set"};
    }
    for (auto const& name : settings.papi_events) {
        auto event = PAPI_NULL;
        if (PAPI_event_name_to_code(name.c_str(), &event) != PAPI_OK) {
            throw std::runtime_error{"Failed to resolve PAPI event '" + name + '\''};
        }
        if (PAPI_add_event(event_set, event) != PAPI_OK) {
            throw std::runtime_error{"Failed to add PAPI event '" + name + '\''};
        }
    }
    return event_set;
}
#endif

void benchmark_thread(Context context) {
#ifdef WITH_PAPI
    if (!context.settings().papi_events.empty()) {
        context.thread_data().push_papi_event_counter.resize(context.settings().papi_events.size());
        context.thread_data().pop_papi_event_counter.resize(context.settings().papi_events.size());
        try {
            context.thread_data().event_set = prepare_papi(context.settings());
        } catch (std::exception const& e) {
            context.write(std::cerr) << e.what() << '\n';
        }
    }
#endif
    std::vector<key_type> prefill(static_cast<std::size_t>(context.settings().prefill_per_thread));
    if (context.id() == 0) {
        std::clog << "Preparing...\n";
    }
    std::seed_seq seed{context.settings().seed, context.id()};
    std::default_random_engine rng(seed);
    context.synchronize();
    std::generate(
        prefill.begin(), prefill.end(),
        [&rng, min = 1UL,
         max = static_cast<key_type>(context.settings().elements_per_thread * context.settings().num_threads)]() {
            return std::uniform_int_distribution<key_type>(min, max)(rng);
        });
    std::generate_n(
        context.shared_data().keys.begin() + context.id() * context.settings().elements_per_thread,
        context.settings().elements_per_thread,
        [&rng, min = 1UL,
         max = static_cast<key_type>(context.settings().elements_per_thread * context.settings().num_threads)]() {
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
    SharedData shared_data;
    shared_data.keys.resize(static_cast<std::size_t>(settings.num_threads * settings.elements_per_thread));
    shared_data.thread_data.resize(static_cast<std::size_t>(settings.num_threads));

    auto pq = pq_type(
        settings.num_threads,
        static_cast<std::size_t>(settings.num_threads * (settings.elements_per_thread + settings.prefill_per_thread)),
        settings.pq_settings);

    thread_coordination::dispatch(settings.affinity, settings.num_threads, [&](auto ctx) {
        benchmark_thread(Context(std::move(ctx), pq.get_handle(), shared_data, settings));
    });

    std::clog << "Done\n";
    std::clog << '\n';
    std::clog << "= Results =\n";
    std::clog << "Push time: " << std::fixed << std::setprecision(3)
              << std::chrono::duration<double>(shared_data.push_time).count() << " s\n";
    std::clog << "Pop time: " << std::fixed << std::setprecision(3)
              << std::chrono::duration<double>(shared_data.pop_time).count() << " s\n";
    write_result_json(settings, shared_data, std::cout);
}

int main(int argc, char* argv[]) {
    benchmark::write_run_header<pq_type>(argc, argv, std::clog);

    cxxopts::Options cmd(argv[0]);
    cmd.add_options()("h,help", "Print this help", cxxopts::value<bool>());
    Settings settings{};
    register_cmd_options(settings, cmd);

    try {
        auto args = cmd.parse(argc, argv);
        if (args.count("help") > 0) {
            std::clog << cmd.help() << '\n';
            return EXIT_SUCCESS;
        }
    } catch (std::exception const& e) {
        std::cerr << "Error parsing command line: " << e.what() << '\n';
        std::cerr << "Use --help for usage information" << '\n';
        return EXIT_FAILURE;
    }

    std::clog << "= Settings =\n";
    write_settings_human_readable(settings, std::clog);
    std::clog << '\n';

    if (!validate_settings(settings)) {
        return EXIT_FAILURE;
    }

    std::clog << "= Running benchmark =\n";
    run_benchmark(settings);
    return EXIT_SUCCESS;
}
