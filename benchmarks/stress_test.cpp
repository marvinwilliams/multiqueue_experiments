#include "util/build_info.hpp"
#include "util/task.hpp"
#include "wrapper/selector.hpp"

#include "cxxopts.hpp"

#ifdef LOG_OPERATIONS
#include "util/operation_log.hpp"
#endif

#ifdef WITH_PAPI
#include <papi.h>
#include <pthread.h>
#endif

#include <atomic>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <new>
#include <numeric>
#include <random>
#include <type_traits>
#include <utility>
#include <vector>

using key_type = unsigned long;
using value_type = unsigned long;

using pq_type = PQ<true, key_type, std::pair<key_type, value_type>>;
using handle_type = pq_type::handle_type;

struct Settings {
    enum class Mode { Update, Random, Pop, PushRandom, PushAscending };

    int num_threads = 4;
    long long prefill_per_thread = 1 << 20;
    long long num_iterations = 1 << 24;
    Mode mode = Mode::Update;
    key_type min_prefill = 1;
    key_type max_prefill = 1 << 20;
    long min_update = 0;
    long max_update = 1 << 20;
    key_type min_random = 1;
    key_type max_random = 1 << 20;
    long long batch_size = 1 << 12;
    int seed = 1;
    std::chrono::seconds timeout{0};
#ifdef LOG_OPERATIONS
    std::filesystem::path log_file{"operation_log.txt"};
#endif
#ifdef WITH_PAPI
    std::vector<std::string> papi_events;
#endif

    static constexpr int pushes_per_iteration(Mode mode) noexcept {
        return mode == Mode::Pop ? 0 : 1;
    }

    static constexpr int pops_per_iteration(Mode mode) noexcept {
        return mode == Mode::PushRandom || mode == Mode::PushAscending ? 0 : 1;
    }

    static bool validate(Settings const& settings) {
        if (settings.num_threads <= 0) {
            std::cerr << "Error: Number of threads must be greater than 0\n";
            return false;
        }
        if (settings.min_prefill <= 0 || settings.min_random <= 0) {
            std::cerr << "Error: Keys must be greater than 0\n";
            return false;
        }
        if (settings.max_prefill < settings.min_prefill || settings.max_random < settings.min_random ||
            settings.max_update < settings.min_update) {
            std::cerr << "Error: Max must be greater than min\n";
            return false;
        }
        if (settings.batch_size <= 0) {
            std::cerr << "Error: batch size must be greater than 0\n";
            return false;
        }
        if (settings.prefill_per_thread < 0 || settings.num_iterations < 0) {
            std::cerr << "Error: Prefill and iterations must be nonnegative\n";
            return false;
        }
        if (settings.num_iterations * pops_per_iteration(settings.mode) >
            settings.num_threads * settings.prefill_per_thread +
                settings.num_iterations * pushes_per_iteration(settings.mode)) {
            std::cerr << "Error: Number of pops must not be greater than number of pushes\n";
            return false;
        }
#ifdef WITH_PAPI
        for (auto const& name : settings.papi_events) {
            if (PAPI_query_named_event(name.c_str()) != PAPI_OK) {
                std::cerr << "Error: PAPI event " << name << " not available\n";
                return false;
            }
        }
#endif
        return true;
    }
};

struct ThreadData {
    std::chrono::high_resolution_clock::time_point start_time{};
    std::chrono::high_resolution_clock::time_point end_time{};
    long long num_iterations = 0;
    long long num_failed_pops = 0;
#ifdef WITH_PAPI
    bool papi_success = true;
    std::vector<long long> papi_counter{};
#endif
#ifdef LOG_OPERATIONS
    struct PushLog {
        long long tick;
        unsigned long key;
        std::size_t index;
    };
    struct PopLog {
        long long tick;
        std::size_t ref_index;
    };
    std::vector<PushLog> pushes;
    std::vector<PopLog> pops;
#endif
};

Settings::Mode parse_mode(char c) {
    switch (c) {
        case 'u':
            return Settings::Mode::Update;
        case 'r':
            return Settings::Mode::Random;
        case 'p':
            return Settings::Mode::Pop;
        case 's':
            return Settings::Mode::PushRandom;
        case 'a':
            return Settings::Mode::PushAscending;
        default:
            throw std::invalid_argument("Invalid work mode");
    }
}

static constexpr auto mode_name(Settings::Mode mode) noexcept {
    switch (mode) {
        case Settings::Mode::Update:
            return "update";
        case Settings::Mode::Random:
            return "random";
        case Settings::Mode::Pop:
            return "pop";
        case Settings::Mode::PushRandom:
            return "push";
        case Settings::Mode::PushAscending:
            return "push (ascending)";
        default:
            return "unknown";
    }
};

void write_settings(Settings const& settings, std::ostream& out) {
    out << "Threads: " << settings.num_threads << '\n'
        << "Prefill per thread: " << settings.prefill_per_thread << '\n'
        << "Prefill range: [" << static_cast<unsigned long>(settings.min_prefill) << ", "
        << static_cast<unsigned long>(settings.max_prefill) << "]\n"
        << "Iterations: " << settings.num_iterations << '\n'
        << "Mode: " << mode_name(settings.mode) << '\n';
    if (settings.mode == Settings::Mode::Update) {
        out << "Update range: [" << settings.min_update << ", " << settings.max_update << "]\n";
    } else if (settings.mode == Settings::Mode::Random || settings.mode == Settings::Mode::PushRandom) {
        out << "Key range: [" << static_cast<unsigned long>(settings.min_random) << ", "
            << static_cast<unsigned long>(settings.max_random) << "]\n";
    }
    out << "batch size: " << settings.batch_size << '\n' << "Seed: " << settings.seed << '\n';
    if (settings.timeout.count() > 0) {
        out << "Timeout (ms): " << settings.timeout.count() << "ms\n";
    }
#ifdef WITH_PAPI
    if (!settings.papi_events.empty()) {
        out << "PAPI events: ";
        std::copy(settings.papi_events.begin(), settings.papi_events.end(),
                  std::ostream_iterator<std::string>(out, " "));
        out << '\n';
    }
#endif
#ifdef LOG_OPERATIONS
    out << "Log operations to: " << settings.log_file << '\n';
#endif
}

void write_result_json(Settings const& settings, std::vector<ThreadData> const& data, std::ostream& out) {
    auto write_vector = [&out](auto const& v, auto&& f) {
        out << '[';
        for (auto it = v.begin(); it != std::prev(v.end()); ++it) {
            out << f(*it) << ',';
        }
        out << f(*std::prev(v.end())) << ']';
    };
    auto last_end = std::max_element(data.begin(), data.end(), [](auto const& lhs, auto const& rhs) {
                        return lhs.end_time < rhs.end_time;
                    })->end_time;
    auto first_start = std::min_element(data.begin(), data.end(), [](auto const& lhs, auto const& rhs) {
                           return lhs.start_time < rhs.start_time;
                       })->start_time;
    out << '{' << std::quoted("threads") << ':' << settings.num_threads << ',' << std::quoted("prefill") << ':'
        << settings.prefill_per_thread << ',' << std::quoted("prefill_min") << ':' << settings.min_prefill << ','
        << std::quoted("prefill_max") << ':' << settings.max_prefill << ',' << std::quoted("iterations") << ':'
        << settings.num_iterations << ',' << std::quoted("mode") << ':' << std::quoted(mode_name(settings.mode));
    if (settings.mode == Settings::Mode::Update) {
        out << ',' << std::quoted("update_min") << ':' << settings.min_update << ',' << std::quoted("update_max") << ':'
            << settings.max_update;
    } else if (settings.mode == Settings::Mode::Random || settings.mode == Settings::Mode::PushRandom) {
        out << ',' << std::quoted("random_min") << ':' << settings.min_random << ',' << std::quoted("random_max") << ':'
            << settings.max_random;
    }
    out << ',' << std::quoted("batch_size") << ':' << settings.batch_size << ',' << std::quoted("seed") << ':'
        << settings.seed << ',' << std::quoted("timeout") << ':' << settings.timeout.count() << ','
        << std::quoted("total_time") << ':' << (last_end - first_start).count() << ',' << std::quoted("thread_time")
        << ':';
    write_vector(data, [](auto const& e) { return (e.end_time - e.start_time).count(); });
    out << ',' << std::quoted("thread_pushes") << ':';
    write_vector(data,
                 [i = Settings::pushes_per_iteration(settings.mode)](auto const& e) { return e.num_iterations * i; });
    out << ',' << std::quoted("thread_pops") << ':';
    write_vector(data,
                 [i = Settings::pops_per_iteration(settings.mode)](auto const& e) { return e.num_iterations * i; });
    out << ',' << std::quoted("thread_failed_pops") << ':';
    write_vector(data, [](auto const& e) { return e.num_failed_pops; });
#ifdef WITH_PAPI
    for (std::size_t i = 0; i < settings.papi_events.size(); ++i) {
        out << ',' << settings.papi_events[i] << ':';
        write_vector(data, [i](auto const& e) { return e.papi_counter[i]; });
    }
#endif
    out << '}' << std::endl;
}

void generate_work(Settings const& settings, int id, std::vector<key_type>& prefill, std::vector<long long>& work) {
    std::seed_seq seed{settings.seed, id};
    std::default_random_engine rng(seed);

    prefill.resize(static_cast<std::size_t>(settings.prefill_per_thread));
    std::generate(prefill.begin(), prefill.end(), [&rng, min = settings.min_prefill, max = settings.max_prefill]() {
        return std::uniform_int_distribution<key_type>(min, max)(rng);
    });
    auto work_size_per_thread =
        settings.num_iterations * Settings::pushes_per_iteration(settings.mode) / settings.num_threads;
    auto start = id * work_size_per_thread;
    if (id == settings.num_threads - 1) {
        work_size_per_thread = settings.num_iterations * Settings::pushes_per_iteration(settings.mode) - start;
    }
    switch (settings.mode) {
        case Settings::Mode::Update: {
            std::generate_n(work.begin() + start, work_size_per_thread,
                            [&rng, min = settings.min_update, max = settings.max_update]() {
                                return std::uniform_int_distribution<long>(min, max)(rng);
                            });
            break;
        }
        case Settings::Mode::Random:
        case Settings::Mode::PushRandom: {
            std::generate_n(work.begin() + start, work_size_per_thread,
                            [&rng, min = settings.min_random, max = settings.max_random]() {
                                return static_cast<long>(std::uniform_int_distribution<key_type>(min, max)(rng));
                            });
            break;
        }
        case Settings::Mode::PushAscending: {
            std::generate_n(work.begin() + start, work_size_per_thread, [i = start]() mutable { return i++; });
            break;
        }
        case Settings::Mode::Pop:
            break;
    }
}

#ifdef WITH_PAPI
int start_papi(std::vector<std::string> const& events) {
    int event_set = PAPI_NULL;
    if (int ret = PAPI_register_thread(); ret != PAPI_OK) {
        throw std::runtime_error("Failed to register thread for PAPI");
    }
    if (int ret = PAPI_create_eventset(&event_set); ret != PAPI_OK) {
        throw std::runtime_error("Failed to create PAPI event set");
    }
    for (auto const& name : events) {
        auto event = PAPI_NULL;
        if (PAPI_event_name_to_code(name.c_str(), &event) != PAPI_OK) {
            throw std::runtime_error("Failed to get PAPI event code for event '" + name + '\'');
        }
        if (PAPI_add_event(event_set, event) != PAPI_OK) {
            throw std::runtime_error("Failed to add PAPI event '" + name + '\'');
        }
    }
    if (int ret = PAPI_start(event_set); ret != PAPI_OK) {
        throw std::runtime_error("Failed to start PAPI");
    }
    return event_set;
}
#endif

template <typename Work>
void work_loop(Settings const& settings, task::Control const& tc, ThreadData& data, Work&& work) {
    static std::atomic<long long> counter{0};
    tc.synchronize();
#ifdef WITH_PAPI
    int event_set = PAPI_NULL;
    bool papi_started = false;
    if (!settings.papi_events.empty()) {
        data.papi_counter.resize(settings.papi_events.size());
        try {
            event_set = start_papi(settings.papi_events);
            papi_started = true;
        } catch (std::runtime_error const& e) {
            tc.write(std::cerr) << "Error: " << e.what() << '\n';
            data.papi_success = false;
        }
    }
#endif
    data.start_time = std::chrono::high_resolution_clock::now();
    auto from = counter.fetch_add(settings.batch_size, std::memory_order_relaxed);
    for (; from < settings.num_iterations - settings.batch_size;
         from = counter.fetch_add(settings.batch_size, std::memory_order_relaxed)) {
        for (auto i = 0; i < settings.batch_size; ++i) {
            work(from + i);
        }
        data.num_iterations += settings.batch_size;
        if (settings.timeout > std::chrono::seconds::zero()) {
            auto now = std::chrono::high_resolution_clock::now();
            if (now - data.start_time > settings.timeout) {
                from = settings.num_iterations;
                break;
            }
        }
    }
    for (; from < settings.num_iterations; ++from) {
        work(from);
        ++data.num_iterations;
    }
    data.end_time = std::chrono::high_resolution_clock::now();
#ifdef WITH_PAPI
    if (papi_started) {
        if (int ret = PAPI_stop(event_set, data.papi_counter.data()); ret != PAPI_OK) {
            tc.write(std::cerr) << "Error: Failed to stop performance counters\n";
            data.papi_success = false;
        }
    }
#endif
    tc.synchronize();
}

struct Context {
    std::vector<std::vector<key_type>> prefill;
    std::vector<long long> work;
    std::chrono::steady_clock::time_point start_time;

    std::ostream& log(std::ostream& out = std::clog) const {
        return out << '[' << std::fixed << std::setprecision(2)
                   << std::chrono::duration<double>(std::chrono::steady_clock::now() - start_time).count() << "s] ";
    }
};

ThreadData benchmark_thread(Settings const& settings, task::Control tc, pq_type& pq, Context& context) {
    handle_type handle = pq.get_handle();
    ThreadData data{};
    auto const id = tc.id();
#ifdef LOG_OPERATIONS
    data.pushes.reserve(static_cast<std::size_t>(
        settings.prefill_per_thread + settings.num_iterations * Settings::pushes_per_iteration(settings.mode)));
    data.pops.reserve(static_cast<std::size_t>(settings.num_iterations * Settings::pops_per_iteration(settings.mode)));

    auto pq_push = [&](key_type key, value_type value) {
        handle.push({key, value});
        auto tick = static_cast<long long>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        data.pushes.push_back({tick, key, static_cast<std::size_t>(value)});
    };

    auto pq_pop = [&]() {
        while (true) {
            auto tick = static_cast<long long>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
            auto retval = handle.try_pop();
            if (retval) {
                data.pops.push_back({tick, static_cast<std::size_t>(retval->second)});
                return retval->first;
            }
        }
    };
#else
    auto pq_push = [&](key_type key, value_type value) { handle.push({key, value}); };
    auto pq_pop = [&]() {
        while (true) {
            auto retval = handle.try_pop();
            if (retval) {
                return retval->first;
            }
        }
    };
#endif

    tc.once([&]() { context.log() << "Generating work...\n"; });
    tc.synchronize();
    generate_work(settings, id, context.prefill[static_cast<std::size_t>(id)], context.work);
    tc.synchronize();
    tc.once([&]() { context.log() << "Prefilling...\n"; });
    tc.synchronize();
    for (auto i = 0LL; i < settings.prefill_per_thread; ++i) {
        pq_push(context.prefill[static_cast<std::size_t>(id)][static_cast<std::size_t>(i)],
                static_cast<value_type>(id * settings.prefill_per_thread + i));
    }
    tc.synchronize();
    tc.once([&]() { context.log() << "Running benchmark...\n"; });
    switch (settings.mode) {
        case Settings::Mode::Update: {
            work_loop(settings, tc, data, [&](int i) {
                auto key = pq_pop();
                auto new_key =
                    static_cast<key_type>(static_cast<long long>(key) + context.work[static_cast<std::size_t>(i)]);
                pq_push(new_key, static_cast<value_type>(settings.num_threads * settings.prefill_per_thread + i));
            });
            break;
        }
        case Settings::Mode::Random: {
            work_loop(settings, tc, data, [&](int i) {
                pq_pop();
                auto new_key = static_cast<key_type>(context.work[static_cast<std::size_t>(i)]);
                pq_push(new_key, static_cast<value_type>(settings.num_threads * settings.prefill_per_thread + i));
            });
            break;
        }
        case Settings::Mode::PushRandom:
        case Settings::Mode::PushAscending: {
            work_loop(settings, tc, data, [&](int i) {
                auto new_key = static_cast<key_type>(context.work[static_cast<std::size_t>(i)]);
                pq_push(new_key, static_cast<value_type>(settings.num_threads * settings.prefill_per_thread + i));
            });
            break;
        }
        case Settings::Mode::Pop: {
            work_loop(settings, tc, data, [&](int) { pq_pop(); });
            break;
        }
    }
    return data;
}

bool run_benchmark(Settings const& settings, cxxopts::ParseResult const& parse_result) {
#ifdef LOG_OPERATIONS
    std::ofstream log_out;
    if (!settings.log_file.empty()) {
        log_out = std::ofstream(settings.log_file);
    }
    if (!log_out) {
        std::cerr << "Error: Could not open file " << settings.log_file << " for writing" << std::endl;
        return false;
    }
#endif

    Context context{};
    context.start_time = std::chrono::steady_clock::now();

    context.prefill.resize(static_cast<std::size_t>(settings.num_threads));
    context.work.resize(
        static_cast<std::size_t>(settings.num_iterations * Settings::pushes_per_iteration(settings.mode)));
    auto max_capacity = context.prefill.size() +
        (settings.mode == Settings::Mode::PushAscending || settings.mode == Settings::Mode::PushRandom
             ? context.work.size()
             : 0);
    auto pq = pq_type(settings.num_threads, 2 * max_capacity, parse_result);
    std::clog << "Priority queue: ";
    pq.describe(std::clog) << '\n' << '\n';
    std::vector<ThreadData> thread_data(static_cast<std::size_t>(settings.num_threads));
    task::Runner runner(affinity::NUMA{cores_per_numa_node, num_numa_nodes}, settings.num_threads, [&](auto tc) {
        thread_data[static_cast<std::size_t>(tc.id())] = benchmark_thread(settings, tc, pq, context);
    });
    runner.wait();

#ifdef LOG_OPERATIONS
    std::vector<ThreadData::PushLog> pushes;
    pushes.reserve(std::accumulate(thread_data.begin(), thread_data.end(), 0UL,
                                   [](std::size_t sum, auto const& e) { return sum + e.pushes.size(); }));
    std::vector<ThreadData::PopLog> pops;
    pushes.reserve(std::accumulate(thread_data.begin(), thread_data.end(), 0UL,
                                   [](std::size_t sum, auto const& e) { return sum + e.pops.size(); }));
    for (auto const& e : thread_data) {
        pushes.insert(pushes.end(), e.pushes.begin(), e.pushes.end());
        pops.insert(pops.end(), e.pops.begin(), e.pops.end());
    }
    std::sort(pushes.begin(), pushes.end(), [](auto const& lhs, auto const& rhs) { return lhs.tick < rhs.tick; });
    std::vector<std::size_t> push_index(pushes.size());
    for (std::size_t i = 0; i < pushes.size(); ++i) {
        push_index[pushes[i].index] = i;
    }
    std::sort(pops.begin(), pops.end(), [](auto const& lhs, auto const& rhs) { return lhs.tick < rhs.tick; });
    auto push_it = pushes.begin();
    for (auto const& pop : pops) {
        auto element_it = std::next(pushes.begin(), static_cast<std::ptrdiff_t>(push_index[pop.ref_index]));
        while (push_it != pushes.end() && (push_it->tick < pop.tick || push_it <= element_it)) {
            log_out << push_it->key << '\n';
            ++push_it;
        }
        log_out << -static_cast<long long>(push_index[pop.ref_index]) << '\n';
    }
    while (push_it != pushes.end()) {
        log_out << push_it->key << '\n';
        ++push_it;
    }
#endif
    context.log() << "Finished\n";
    write_result_json(settings, thread_data, std::cout);
    return true;
}

int main(int argc, char* argv[]) {
    write_build_info(std::clog);
    Settings settings;
    cxxopts::Options cmd(argv[0]);
    int timeout_ms = 0;
    cmd.add_options()
        // clang-format off
        ("h,help", "Print this help", cxxopts::value<bool>())
        ("j,threads", "Number of threads", cxxopts::value<int>(settings.num_threads), "NUMBER")
        ("p,prefill", "Prefill per thread", cxxopts::value<long long>(settings.prefill_per_thread), "NUMBER")
        ("n,iterations", "Number of iterations", cxxopts::value<long long>(settings.num_iterations), "NUMBER")
        ("mode", "Operation mode ([u]pdate, [r]andom, [p]op, pu[s]h, push [a]scending", cxxopts::value<char>(), "STRING")
        ("min-prefill", "Min prefill key", cxxopts::value<key_type>(settings.min_prefill), "NUMBER")
        ("max-prefill", "Max prefill key", cxxopts::value<key_type>(settings.max_prefill), "NUMBER")
        ("min-update", "Min update", cxxopts::value<long>(settings.min_update), "NUMBER")
        ("max-update", "Max update", cxxopts::value<long>(settings.max_update), "NUMBER")
        ("min-random", "Min random key", cxxopts::value<key_type>(settings.min_random), "NUMBER")
        ("max-random", "Max random key", cxxopts::value<key_type>(settings.max_random), "NUMBER")
        ("batch-size", "Batch size", cxxopts::value<long long>(settings.batch_size), "NUMBER")
        ("s,seed", "Initial seed", cxxopts::value<int>(settings.seed), "NUMBER")
        ("timeout", "Timeout in milliseconds", cxxopts::value<int>(timeout_ms), "NUMBER")
#ifdef LOG_OPERATIONS
        ("log-file", "File to write the operation log to", cxxopts::value<std::filesystem::path>(settings.log_file), "PATH")
#endif
#ifdef WITH_PAPI
        ("r,pc", "Performance counters", cxxopts::value<std::vector<std::string>>(settings.papi_events))
#endif
        // clang-format on
        ;
    pq_type::add_options(cmd);
#ifdef LOG_OPERATIONS
    cmd.parse_positional({"log-file"});
#endif
    cxxopts::ParseResult args;
    try {
        args = cmd.parse(argc, argv);
    } catch (std::exception const& e) {
        std::cerr << "Error parsing command line: " << e.what() << std::endl;
        std::cerr << cmd.help() << std::endl;
        return EXIT_FAILURE;
    }

    if (args.count("help") > 0) {
        std::clog << cmd.help() << std::endl;
        return EXIT_SUCCESS;
    }
    if (args.count("mode") > 0) {
        try {
            settings.mode = parse_mode(args["mode"].as<char>());
        } catch (std::invalid_argument const& e) {
            std::cerr << "Error parsing command line: " << e.what() << std::endl;
            std::cerr << cmd.help() << std::endl;
            return EXIT_FAILURE;
        }
    }
    settings.timeout = std::chrono::seconds(timeout_ms);

    std::clog << '\n';
    std::clog << "Command line:";
    for (int i = 0; i < argc; ++i) {
        std::clog << ' ' << argv[i];
    }
    std::clog << '\n';
    write_settings(settings, std::clog);

#ifdef WITH_PAPI
    if (!settings.papi_events.empty()) {
        if (int ret = PAPI_library_init(PAPI_VER_CURRENT); ret != PAPI_VER_CURRENT) {
            std::cerr << "Error: Failed to initialize PAPI library\n";
            return EXIT_FAILURE;
        }
        if (int ret = PAPI_thread_init(pthread_self); ret != PAPI_OK) {
            std::cerr << "Error: Failed to initialize PAPI thread support\n";
            return EXIT_FAILURE;
        }
    }
#endif

    if (!Settings::validate(settings)) {
        return EXIT_FAILURE;
    }

    bool success = run_benchmark(settings, args);
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
