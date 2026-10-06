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
#include <type_traits>
#include <utility>
#include <vector>

#ifdef __SSE2__
#include <emmintrin.h>
#define PAUSE _mm_pause()
#else
#define PAUSE void(0)
#endif

using key_type = unsigned long;
using value_type = unsigned long;

using pq_type = PQ<true, key_type, value_type>;
using handle_type = pq_type::handle_type;

struct Settings {
    benchmark::BaseSettings base_settings{};
    pq_type::settings_type pq_settings;
    long long prefill_per_thread = 1 << 20;
    long long iterations_per_thread = 1 << 24;
    key_type min_prefill = 1;
    key_type max_prefill = 1 << 20;
    long min_update = 0;
    long max_update = 1 << 20;
    long long batch_size = 1 << 12;
    int seed = 1;
    int timeout_s = 0;
    int sleep_us = 0;
#ifdef LOG_OPERATIONS
    std::filesystem::path log_file = "operation_log.txt";
#endif
#ifdef WITH_PAPI
    std::vector<std::string> papi_events;
#endif

    void register_cmd_options(cxxopts::Options& cmd) {
        base_settings.register_cmd_options(cmd);
        pq_settings.register_cmd_options(cmd);
        cmd.add_options()
            // clang-format off
            ("p,prefill", "Prefill per thread", cxxopts::value<long long>(prefill_per_thread), "NUMBER")
            ("n,iterations", "Number of iterations per thread", cxxopts::value<long long>(iterations_per_thread), "NUMBER")
            ("min-prefill", "Min prefill key", cxxopts::value<key_type>(min_prefill), "NUMBER")
            ("max-prefill", "Max prefill key", cxxopts::value<key_type>(max_prefill), "NUMBER")
            ("min-update", "Min update", cxxopts::value<long>(min_update), "NUMBER")
            ("max-update", "Max update", cxxopts::value<long>(max_update), "NUMBER")
            ("batch-size", "Batch size", cxxopts::value<long long>(batch_size), "NUMBER")
            ("s,seed", "Initial seed", cxxopts::value<int>(seed), "NUMBER")
            ("t,timeout", "Timeout in seconds", cxxopts::value<int>(timeout_s), "NUMBER")
            ("q,sleep", "Time in microseconds to wait between operations", cxxopts::value<int>(sleep_us), "NUMBER")
#ifdef LOG_OPERATIONS
            ("l,log-file", "File to write the operation log to", cxxopts::value<std::filesystem::path>(log_file), "PATH")
#endif
#ifdef WITH_PAPI
            ("r,count-event", "Papi event to count", cxxopts::value<std::vector<std::string>>(papi_events), "STRING")
#endif
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
        if (iterations_per_thread < 0) {
            std::cerr << "Error: Iterations must be nonnegative\n";
            return false;
        }
        if (min_prefill <= 0) {
            std::cerr << "Error: Prefill keys must be greater than 0\n";
            return false;
        }
        if (max_prefill < min_prefill) {
            std::cerr << "Error: Invalid prefill range\n";
            return false;
        }
        if (min_update < 0) {
            std::cerr << "Error: Min update must be nonnegative\n";
            return false;
        }
        if (max_update < min_update) {
            std::cerr << "Error: Invalid update range\n";
            return false;
        }
        if (batch_size <= 0) {
            std::cerr << "Error: batch size must be greater than 0\n";
            return false;
        }
        if (timeout_s < 0) {
            std::cerr << "Error: Timeout must be nonnegative\n";
            return false;
        }
        if (sleep_us < 0) {
            std::cerr << "Error: Sleep must be nonnegative\n";
            return false;
        }
        if (seed <= 0) {
            std::cerr << "Error: Seed must be greater than 0\n";
            return false;
        }
#ifdef LOG_OPERATIONS
        if (log_file.empty()) {
            std::cerr << "Error: Log file name must not be empty\n";
            return false;
        }
        auto out = std::ofstream(log_file);
        if (out.fail()) {
            std::cerr << "Error: Could not open file " << log_file << " for writing\n";
            return false;
        }
        out.close();
#endif
#ifdef WITH_PAPI
        if (!papi_events.empty()) {
            if (int ret = PAPI_library_init(PAPI_VER_CURRENT); ret != PAPI_VER_CURRENT) {
                std::cerr << "Error: Failed to initialize PAPI library\n";
                return false;
            }
            if (int ret = PAPI_thread_init(pthread_self); ret != PAPI_OK) {
                std::cerr << "Error: Failed to initialize PAPI thread support\n";
                return false;
            }
            for (auto const& name : papi_events) {
                if (PAPI_query_named_event(name.c_str()) != PAPI_OK) {
                    std::cerr << "Error: PAPI event '" << name << "' not available\n";
                    return false;
                }
            }
        }
#endif
        return true;
    }

    void write_human_readable(std::ostream& out) const {
        base_settings.write_human_readable(out);
        pq_settings.write_human_readable(out);
        out << "Prefill per thread: " << prefill_per_thread << '\n';
        out << "Iterations per thread: " << iterations_per_thread << '\n';
        out << "Prefill range: [" << min_prefill << ", " << max_prefill << "]\n";
        out << "Update range: [" << min_update << ", " << max_update << "]\n";
        out << "Batch size: " << batch_size << '\n';
        out << "Timeout: ";
        if (timeout_s == 0) {
            out << "None\n";
        } else {
            out << timeout_s << " s\n";
        }
        out << "Sleep: ";
        if (sleep_us == 0) {
            out << "None\n";
        } else {
            out << sleep_us << " us\n";
        }
        out << "Seed: " << seed << '\n';
#ifdef LOG_OPERATIONS
        out << "Log file: " << log_file << '\n';
#endif
#ifdef WITH_PAPI
        out << "PAPI events: [";
        for (std::size_t i = 0; i < papi_events.size(); ++i) {
            out << papi_events[i];
            if (i != papi_events.size() - 1) {
                out << ", ";
            }
        }
        out << ']' << '\n';
#endif
    }

    void write_json(json::Object& obj) const {
        base_settings.write_json(obj);
        obj.object("pq", [this](json::Object& o) { pq_settings.write_json(o); });
        obj.entry("prefill_per_thread", prefill_per_thread);
        obj.entry("iterations_per_thread", iterations_per_thread);
        obj.entry("prefill_min", min_prefill);
        obj.entry("prefill_max", max_prefill);
        obj.entry("update_min", min_update);
        obj.entry("update_max", max_update);
        obj.entry("batch_size", batch_size);
        obj.entry("timeout_s", timeout_s);
        obj.entry("sleep_us", sleep_us);
        obj.entry("seed", seed);
#ifdef WITH_PAPI
        obj.array("papi_events", papi_events);
#endif
    }
};

struct ThreadData {
    long long iter_count = 0;
    long long failed_pop_count = 0;
#ifdef WITH_PAPI
    std::vector<long long> papi_event_counter{};
#endif
#ifdef LOG_OPERATIONS
    struct PushLog {
        std::chrono::high_resolution_clock::time_point tick;
        std::pair<key_type, value_type> element;
    };
    struct PopLog {
        std::chrono::high_resolution_clock::time_point tick;
        value_type val;
    };
    std::vector<PushLog> pushes;
    std::vector<PopLog> pops;
#endif

    void write_json(json::Object& obj) const {
        obj.entry("iterations", iter_count);
        obj.entry("failed_pops", failed_pop_count);
#ifdef WITH_PAPI
        obj.array("papi_event_counter", papi_event_counter);
#endif
    }
};

#ifdef LOG_OPERATIONS
void write_log(std::vector<ThreadData> const& thread_data, std::ostream& out) {
    std::vector<ThreadData::PushLog> pushes;
    pushes.reserve(std::accumulate(thread_data.begin(), thread_data.end(), 0UL,
                                   [](std::size_t sum, auto const& e) { return sum + e.pushes.size(); }));
    std::vector<ThreadData::PopLog> pops;
    pops.reserve(std::accumulate(thread_data.begin(), thread_data.end(), 0UL,
                                 [](std::size_t sum, auto const& e) { return sum + e.pops.size(); }));
    for (auto const& e : thread_data) {
        pushes.insert(pushes.end(), e.pushes.begin(), e.pushes.end());
        pops.insert(pops.end(), e.pops.begin(), e.pops.end());
    }
    std::sort(pushes.begin(), pushes.end(), [](auto const& lhs, auto const& rhs) { return lhs.tick < rhs.tick; });
    std::vector<std::size_t> push_index(pushes.size());
    for (std::size_t i = 0; i < pushes.size(); ++i) {
        push_index[pushes[i].element.second] = i;
    }
    std::sort(pops.begin(), pops.end(), [](auto const& lhs, auto const& rhs) { return lhs.tick < rhs.tick; });
    out << pushes.size() << ' ' << pops.size() << '\n';
    std::size_t i = 0;
    for (auto const& pop : pops) {
        while ((i != pushes.size() && pushes[i].tick < pop.tick)) {
            out << '+' << pushes[i].element.first << '\n';
            ++i;
        }
        out << '-' << push_index[static_cast<std::size_t>(pop.val)] << '\n';
    }
    for (; i < pushes.size(); ++i) {
        out << '+' << pushes[i].element.first << '\n';
    }
}
#endif

struct SharedData {
    std::vector<long long> updates;
    std::atomic_llong counter{0};
    std::chrono::steady_clock::time_point start_time;
    std::chrono::steady_clock::time_point end_time;
    std::vector<ThreadData> thread_data;
};

void write_result_json(Settings const& settings, SharedData const& data, std::ostream& out) {
    {
        json::Object root{out};
        root.object("settings", [&settings](json::Object& obj) { settings.write_json(obj); });
        root.object("results", [&data](json::Object& results) {
            results.entry("time_ns", std::chrono::nanoseconds{data.end_time - data.start_time}.count());
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

#ifdef LOG_OPERATIONS
    void push(std::pair<key_type, value_type> const& e) {
        handle_.push(e);
        auto tick = std::chrono::high_resolution_clock::now();
        thread_data_.pushes.push_back({tick, e});
    }

    auto try_pop() {
        auto tick = std::chrono::high_resolution_clock::now();
        auto retval = handle_.try_pop();
        if (retval) {
            thread_data_.pops.push_back({tick, retval->second});
        }
        return retval;
    }
#else
    void push(std::pair<key_type, value_type> const& e) {
        handle_.push(e);
    }

    auto try_pop() {
        return handle_.try_pop();
    }
#endif

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

[[gnu::noinline]] void work_loop(Context& context) {
    auto timeout = [t = std::chrono::seconds{context.settings().timeout_s},
                    start = context.shared_data().start_time]() {
        if (t == std::chrono::seconds::zero()) {
            return false;
        }
        return std::chrono::steady_clock::now() > start + t;
    };
    auto sleep = [s = std::chrono::microseconds{context.settings().sleep_us}]() {
        if (s == std::chrono::microseconds::zero()) {
            return;
        }
        auto now = std::chrono::steady_clock::now();
        auto sleep_until = now + s;
        do {
            PAUSE;
            now = std::chrono::steady_clock::now();
        } while (now < sleep_until);
    };
    auto offset =
        static_cast<value_type>(context.settings().base_settings.num_threads * context.settings().prefill_per_thread);
    long long max = context.settings().iterations_per_thread * context.settings().base_settings.num_threads;
    for (auto from = context.shared_data().counter.fetch_add(context.settings().batch_size, std::memory_order_relaxed);
         from < max;
         from = context.shared_data().counter.fetch_add(context.settings().batch_size, std::memory_order_relaxed)) {
        auto to = std::min(from + context.settings().batch_size, max);
        for (auto i = from; i < to; ++i) {
            auto e = context.try_pop();
            while (!e) {
                ++context.thread_data().failed_pop_count;
                if (timeout()) {
                    context.thread_data().iter_count += i - from;
                    return;
                }
                e = context.try_pop();
            }
            sleep();
            context.push({static_cast<key_type>(static_cast<long long>(e->first) +
                                                context.shared_data().updates[static_cast<std::size_t>(i)]),
                          offset + static_cast<value_type>(i)});
        }
        context.thread_data().iter_count += to - from;
        if (timeout()) {
            break;
        }
    }
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
    int event_set = PAPI_NULL;
    if (!context.settings().papi_events.empty()) {
        context.thread_data().papi_event_counter.resize(context.settings().papi_events.size());
        try {
            event_set = prepare_papi(context.settings());
        } catch (std::exception const& e) {
            context.write(std::cerr) << e.what() << '\n';
        }
    }
#endif
#ifdef LOG_OPERATIONS
    context.thread_data().pushes.reserve(
        static_cast<std::size_t>(context.settings().prefill_per_thread + 2 * context.settings().iterations_per_thread));
    context.thread_data().pops.reserve(static_cast<std::size_t>(2 * context.settings().iterations_per_thread));
#endif

    std::vector<key_type> prefill(static_cast<std::size_t>(context.settings().prefill_per_thread));

    if (context.id() == 0) {
        std::clog << "Preparing...\n";
    }
    std::seed_seq seed{context.settings().seed, context.id()};
    std::default_random_engine rng(seed);
    context.synchronize();
    std::generate(prefill.begin(), prefill.end(),
                  [&rng, min = context.settings().min_prefill, max = context.settings().max_prefill]() {
                      return std::uniform_int_distribution<key_type>(min, max)(rng);
                  });
    std::generate_n(context.shared_data().updates.begin() + context.id() * context.settings().iterations_per_thread,
                    context.settings().iterations_per_thread,
                    [&rng, min = context.settings().min_update, max = context.settings().max_update]() {
                        return std::uniform_int_distribution<long>(min, max)(rng);
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
        std::clog << "Working...\n";
    }
#ifdef WITH_PAPI
    if (!context.settings().papi_events.empty()) {
        if (int ret = PAPI_start(event_set); ret != PAPI_OK) {
            context.write(std::cerr) << "Failed to start performance counters\n";
        }
    }
#endif
    if (context.id() == 0) {
        context.shared_data().start_time = std::chrono::high_resolution_clock::now();
    }
    context.synchronize();
    work_loop(context);
    context.synchronize();
    if (context.id() == 0) {
        context.shared_data().end_time = std::chrono::high_resolution_clock::now();
    }
#ifdef WITH_PAPI
    if (!context.settings().papi_events.empty()) {
        if (int ret = PAPI_stop(event_set, context.thread_data().papi_event_counter.data()); ret != PAPI_OK) {
            context.write(std::cerr) << "Failed to stop performance counters\n";
        }
    }
#endif
    context.shared_data().thread_data[static_cast<std::size_t>(context.id())] = std::move(context.thread_data());
}

void run_benchmark(Settings const& settings) {
    SharedData shared_data;
    shared_data.updates.resize(static_cast<std::size_t>(settings.iterations_per_thread * settings.num_threads));
    shared_data.thread_data.resize(static_cast<std::size_t>(settings.num_threads));
    auto pq =
        pq_type(settings.num_threads, static_cast<std::size_t>(settings.prefill_per_thread * settings.num_threads),
                settings.pq_settings);

    thread_coordination::dispatch(settings.base_settings.affinity, settings.base_settings.num_threads, [&](auto ctx) {
        benchmark_thread(Context(std::move(ctx), pq.get_handle(), shared_data, settings));
    });

#ifdef LOG_OPERATIONS
    std::clog << "Writing logs...\n";
    std::ofstream log_out(settings.log_file);  // assumed to be valid
    write_log(shared_data.thread_data, log_out);
    log_out.close();
#endif
    std::clog << "Done\n";
    std::clog << '\n';
    std::clog << "= Results =\n";
    std::clog << "Time (s): " << std::fixed << std::setprecision(3)
              << std::chrono::duration<double>(shared_data.end_time - shared_data.start_time).count() << '\n';
    write_result_json(settings, shared_data, std::cout);
}

int main(int argc, char* argv[]) {
    benchmark::write_header<pq_type>(argc, argv, std::clog);

    cxxopts::Options cmd(argv[0]);
    cmd.add_options()("h,help", "Print this help", cxxopts::value<bool>());
    Settings settings{};
    settings.register_cmd_options(cmd);

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
