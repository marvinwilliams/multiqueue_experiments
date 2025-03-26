#include "util/build_info.hpp"
#include "util/selector.hpp"
#include "util/thread_coordination.hpp"
#ifdef LOG_OPERATIONS
#include "util/operation_log.hpp"
#endif

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
    int num_threads = 4;
    long long prefill_per_thread = 1 << 20;
    long long iterations_per_thread = 1 << 24;
    key_type min_prefill = 1;
    key_type max_prefill = 1 << 20;
    long min_update = 0;
    long max_update = 1 << 20;
    long long batch_size = 1 << 12;
    int seed = 1;
    int affinity = 6;
    int timeout_s = 0;
    int sleep_us = 0;
#ifdef LOG_OPERATIONS
    std::filesystem::path log_file = "operation_log.txt";
#endif
#ifdef WITH_PAPI
    std::vector<std::string> papi_events;
#endif
    pq_type::settings_type pq_settings;
};

void register_cmd_options(Settings& settings, cxxopts::Options& cmd) {
    cmd.add_options()
        // clang-format off
            ("j,threads", "Number of threads", cxxopts::value<int>(settings.num_threads), "NUMBER")
            ("p,prefill", "Prefill per thread", cxxopts::value<long long>(settings.prefill_per_thread), "NUMBER")
            ("n,iterations", "Number of iterations per thread", cxxopts::value<long long>(settings.iterations_per_thread), "NUMBER")
            ("min-prefill", "Min prefill key", cxxopts::value<key_type>(settings.min_prefill), "NUMBER")
            ("max-prefill", "Max prefill key", cxxopts::value<key_type>(settings.max_prefill), "NUMBER")
            ("min-update", "Min update", cxxopts::value<long>(settings.min_update), "NUMBER")
            ("max-update", "Max update", cxxopts::value<long>(settings.max_update), "NUMBER")
            ("batch-size", "Batch size", cxxopts::value<long long>(settings.batch_size), "NUMBER")
            ("s,seed", "Initial seed", cxxopts::value<int>(settings.seed), "NUMBER")
            ("a,affinity", "CPU affinity ("
                "0: None, "
                "1: Thread Id, "
                "2: Same, "
                "3: Close caches, "
                "4: Far caches, "
                "5: Close L3 Far L1, "
                "6: Far L1 Close L3)"
                , cxxopts::value<int>(settings.affinity), "NUMBER")
            ("t,timeout", "Timeout in seconds", cxxopts::value<int>(settings.timeout_s), "NUMBER")
            ("q,sleep", "Time in microseconds to wait between operations", cxxopts::value<int>(settings.sleep_us), "NUMBER")
#ifdef LOG_OPERATIONS
            ("l,log-file", "File to write the operation log to", cxxopts::value<std::filesystem::path>(settings.log_file), "PATH")
#endif
#ifdef WITH_PAPI
            ("r,count-event", "Papi event to count", cxxopts::value<std::vector<std::string>>(settings.papi_events), "STRING")
#endif
        // clang-format on
        ;
    settings.pq_settings.register_cmd_options(cmd);
}

bool validate_settings(Settings const& settings) {
    if (settings.num_threads <= 0) {
        std::cerr << "Error: Number of threads must be greater than 0\n";
        return false;
    }
    if (settings.prefill_per_thread < 0) {
        std::cerr << "Error: Prefill must be nonnegative\n";
        return false;
    }
    if (settings.iterations_per_thread < 0) {
        std::cerr << "Error: Iterations must be nonnegative\n";
        return false;
    }
    if (settings.min_prefill <= 0) {
        std::cerr << "Error: Prefill keys must be greater than 0\n";
        return false;
    }
    if (settings.max_prefill < settings.min_prefill) {
        std::cerr << "Error: Invalid prefill range\n";
        return false;
    }
    if (settings.min_update < 0) {
        std::cerr << "Error: Min update must be nonnegative\n";
        return false;
    }
    if (settings.max_update < settings.min_update) {
        std::cerr << "Error: Invalid update range\n";
        return false;
    }
    if (settings.batch_size <= 0) {
        std::cerr << "Error: batch size must be greater than 0\n";
        return false;
    }
    if (settings.affinity < 0 || settings.affinity > 6) {
        std::cerr << "Error: Invalid affinity\n";
        return false;
    }
    if (settings.timeout_s < 0) {
        std::cerr << "Error: Timeout must be nonnegative\n";
        return false;
    }
    if (settings.sleep_us < 0) {
        std::cerr << "Error: Sleep must be nonnegative\n";
        return false;
    }
    if (settings.seed <= 0) {
        std::cerr << "Error: Seed must be greater than 0\n";
        return false;
    }
#ifdef LOG_OPERATIONS
    if (settings.log_file.empty()) {
        std::cerr << "Error: Log file name must not be empty\n";
        return false;
    }
    auto out = std::ofstream(settings.log_file);
    if (out.fail()) {
        std::cerr << "Error: Could not open file " << settings.log_file << " for writing\n";
        return false;
    }
    out.close();
#endif
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
    auto affinity_name = [](int a) {
        switch (a) {
            case 0:
                return "None";
            case 1:
                return "Thread Id";
            case 2:
                return "Same";
            case 3:
                return "Close caches";
            case 4:
                return "Far caches";
            case 5:
                return "Close L3 Far L1";
            case 6:
                return "Far L1 Close L3";
            default:
                return "";
        }
    };
    out << "Threads: " << settings.num_threads << '\n';
    out << "Prefill per thread: " << settings.prefill_per_thread << '\n';
    out << "Iterations per thread: " << settings.iterations_per_thread << '\n';
    out << "Prefill range: [" << settings.min_prefill << ", " << settings.max_prefill << "]\n";
    out << "Update range: [" << settings.min_update << ", " << settings.max_update << "]\n";
    out << "Batch size: " << settings.batch_size << '\n';
    out << "Affinity: " << affinity_name(settings.affinity) << '\n';
    out << "Timeout: ";
    if (settings.timeout_s == 0) {
        out << "None\n";
    } else {
        out << settings.timeout_s << " s\n";
    }
    out << "Sleep: ";
    if (settings.sleep_us == 0) {
        out << "None\n";
    } else {
        out << settings.sleep_us << " us\n";
    }
    out << "Seed: " << settings.seed << '\n';
#ifdef LOG_OPERATIONS
    out << "Log file: " << settings.log_file << '\n';
#endif
#ifdef WITH_PAPI
    out << "PAPI events: [";
    for (std::size_t i = 0; i < settings.papi_events.size(); ++i) {
        out << settings.papi_events[i];
        if (i != settings.papi_events.size() - 1) {
            out << ", ";
        }
    }
    out << ']' << '\n';
#endif
    settings.pq_settings.write_human_readable(out);
}

void write_settings_json(Settings const& settings, std::ostream& out) {
    out << '{';
    out << std::quoted("num_threads") << ':' << settings.num_threads << ',';
    out << std::quoted("prefill_per_thread") << ':' << settings.prefill_per_thread << ',';
    out << std::quoted("iterations_per_thread") << ':' << settings.iterations_per_thread << ',';
    out << std::quoted("prefill_min") << ':' << settings.min_prefill << ',';
    out << std::quoted("prefill_max") << ':' << settings.max_prefill << ',';
    out << std::quoted("update_min") << ':' << settings.min_update << ',';
    out << std::quoted("update_max") << ':' << settings.max_update << ',';
    out << std::quoted("batch_size") << ':' << settings.batch_size << ',';
    out << std::quoted("affinity") << ':' << settings.affinity << ',';
    out << std::quoted("timeout_s") << ':' << settings.timeout_s << ',';
    out << std::quoted("sleep_us") << ':' << settings.sleep_us << ',';
    out << std::quoted("seed") << ':' << settings.seed << ',';
#ifdef WITH_PAPI
    out << std::quoted("papi_events") << ':';
    out << '[';
    for (std::size_t i = 0; i < settings.papi_events.size(); ++i) {
        out << std::quoted(settings.papi_events[i]);
        if (i != settings.papi_events.size() - 1) {
            out << ',';
        }
    }
    out << ']' << ',';
#endif
    out << std::quoted("pq") << ':';
    settings.pq_settings.write_json(out);
    out << '}';
}

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
};

void write_thread_data_json(ThreadData const& data, std::ostream& out) {
    out << '{';
    out << std::quoted("iterations") << ':' << data.iter_count << ',';
    out << std::quoted("failed_pops") << ':' << data.failed_pop_count;
#ifdef WITH_PAPI
    out << ',';
    out << std::quoted("papi_event_counter") << ':';
    out << '[';
    for (std::size_t i = 0; i < data.papi_event_counter.size(); ++i) {
        out << data.papi_event_counter[i];
        if (i != data.papi_event_counter.size() - 1) {
            out << ',';
        }
    }
    out << ']';
#endif
    out << '}';
}

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
    std::chrono::high_resolution_clock::time_point start_time;
    std::chrono::high_resolution_clock::time_point end_time;
    std::vector<ThreadData> thread_data;
};

void write_result_json(Settings const& settings, SharedData const& data, std::ostream& out) {
    out << '{';
    out << std::quoted("settings") << ':';
    write_settings_json(settings, out);
    out << ',';
    out << std::quoted("results") << ':';
    out << '{';
    out << std::quoted("time_ns") << ':' << std::chrono::nanoseconds{data.end_time - data.start_time}.count() << ',';
    out << std::quoted("thread_data") << ':';
    out << '[';
    for (auto it = data.thread_data.begin(); it != data.thread_data.end(); ++it) {
        write_thread_data_json(*it, out);
        if (it != std::prev(data.thread_data.end())) {
            out << ',';
        }
    }
    out << ']';
    out << '}';
    out << '}' << '\n';
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
    auto offset = static_cast<value_type>(context.settings().num_threads * context.settings().prefill_per_thread);
    long long max = context.settings().iterations_per_thread * context.settings().num_threads;
    for (auto from = context.shared_data().counter.fetch_add(context.settings().batch_size, std::memory_order_relaxed);
         from < max;
         from = context.shared_data().counter.fetch_add(context.settings().batch_size, std::memory_order_relaxed)) {
        auto to = std::min(from + context.settings().batch_size, max);
        for (auto i = from; i < to; ++i) {
            while (true) {
                if (auto e = context.try_pop(); e) {
                    if (context.settings().sleep_us != 0) {
                        auto sleep_until = std::chrono::high_resolution_clock::now() +
                            std::chrono::microseconds{context.settings().sleep_us};
                        while (std::chrono::high_resolution_clock::now() < sleep_until) {
                            PAUSE;
                        }
                    }
                    context.push({static_cast<key_type>(static_cast<long long>(e->first) +
                                                        context.shared_data().updates[static_cast<std::size_t>(i)]),
                                  offset + static_cast<value_type>(i)});
                    break;
                }
                ++context.thread_data().failed_pop_count;
            }
        }
        context.thread_data().iter_count += to - from;
        if (context.settings().timeout_s != 0 &&
            std::chrono::high_resolution_clock::now() >
                context.shared_data().start_time + std::chrono::seconds{context.settings().timeout_s}) {
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

    auto dispatch = [&](auto const& affinity) {
        auto dispatcher = thread_coordination::Dispatcher(affinity, settings.num_threads, [&](auto ctx) {
            benchmark_thread(Context(std::move(ctx), pq.get_handle(), shared_data, settings));
        });
        dispatcher.wait();
    };
    switch (settings.affinity) {
        case 0:
            dispatch(thread_coordination::affinity::None{});
            break;
        case 1:
            dispatch(thread_coordination::affinity::ThreadId{});
            break;
        case 2:
            dispatch(thread_coordination::affinity::Same{});
            break;
        case 3:
            dispatch(thread_coordination::affinity::CloseCaches{});
            break;
        case 4:
            dispatch(thread_coordination::affinity::FarCaches{});
            break;
        case 5:
            dispatch(thread_coordination::affinity::CloseL3FarL1{});
            break;
        case 6:
            dispatch(thread_coordination::affinity::FarL1CloseL3{});
            break;
    }

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
    write_build_info(std::clog);
    std::clog << '\n';

    std::clog << "= Priority queue =\n";
    pq_type::write_human_readable(std::clog);
    std::clog << '\n';

    std::clog << "= Command line =\n";
    for (int i = 0; i < argc; ++i) {
        std::clog << argv[i];
        if (i != argc - 1) {
            std::clog << ' ';
        }
    }
    std::clog << '\n' << '\n';

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
