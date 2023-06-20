#include "priority_queue_selection.hpp"
#include "thread_coordination.hpp"

#ifdef QUALITY
#include "operation_log.hpp"
#endif

#include "cxxopts.hpp"

#ifdef WITH_PAPI
extern "C" {
#include <papi.h>
#include <pthread.h>
}
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

void print_header() {
    std::clog << "Built on " << __DATE__ << ' ' << __TIME__ << " with\n";
#if defined(__clang__)
    std::clog << "  Clang " << __clang_version__ << '\n';
#elif defined(__GNUC__)
    std::clog << "  GCC " << __VERSION__ << '\n';
#else
    std::clog << "  Unknown compiler\n";
#endif
#ifdef NDEBUG
    std::clog << "  NDEBUG defined\n";
#else
    std::clog << "  NDEBUG not defined\n";
#endif
#if defined WITH_PAPI
    std::clog << "  PAPI " << PAPI_VER_CURRENT << '\n';
#else
    std::clog << "  PAPI unsupported\n";
#endif
}

#ifdef USE_UINT8
#ifdef QUALITY
#error "QUALITY does not work with UINT8"
#endif
using key_type = std::uint8_t;
using value_type = std::uint8_t;
#else
using key_type = unsigned long;
using value_type = unsigned long;
#endif

using pq_type = PriorityQueue<key_type, value_type, true>;
using handle_type = pq_type::handle_type;

struct CommandLineOptions {
    int num_threads = 4;
    std::size_t prefill_per_thread = 1 << 20;
    std::size_t elements_per_thread = 1 << 24;
    char work_mode = 'm';
    int num_push_threads = 1;
    char key_distribution = 'u';
    pq_type::config_type pq_settings{};
#ifdef USE_UINT8
    key_type max_key = std::numeric_limits<key_type>::max();
#else
    key_type max_key = 1 << 30;
#endif
    int seed = 1;
    std::filesystem::path result_file;
#ifdef QUALITY
    std::filesystem::path histogram_file = "histogram.dat";
    std::filesystem::path log_file;
#endif
#ifdef WITH_PAPI
    bool enable_performance_counter = false;
#endif
    bool help = false;
};

enum class WorkMode { Mixed, Split, Increment };
enum class KeyDistribution { Uniform, Ascending, Descending };

[[nodiscard]] auto work_mode_name(WorkMode work_mode) {
    switch (work_mode) {
        case WorkMode::Mixed:
            return "mixed";
        case WorkMode::Split:
            return "split";
        case WorkMode::Increment:
            return "increment";
    }
    return "unknown";
}

[[nodiscard]] auto key_distribution_name(KeyDistribution key_distribution) {
    switch (key_distribution) {
        case KeyDistribution::Uniform:
            return "uniform";
        case KeyDistribution::Ascending:
            return "ascending";
        case KeyDistribution::Descending:
            return "descending";
    }
    return "unknown";
}

struct Settings {
    int num_threads;
    std::size_t prefill_per_thread;
    std::size_t elements_per_thread;
    WorkMode work_mode;
    int num_push_threads;
    KeyDistribution key_distribution;
    key_type max_key;
    int seed;
    std::filesystem::path result_file;
    pq_type::config_type pq_settings;
#ifdef QUALITY
    std::filesystem::path histogram_file;
    std::filesystem::path log_file;
#endif
#ifdef WITH_PAPI
    bool enable_performance_counter;
#endif
    explicit Settings(CommandLineOptions const& command_line_options)
        : num_threads(command_line_options.num_threads),
          prefill_per_thread(command_line_options.prefill_per_thread),
          elements_per_thread(command_line_options.elements_per_thread),
          work_mode(parse_work_mode(command_line_options.work_mode)),
          num_push_threads(command_line_options.num_push_threads),
          key_distribution(parse_key_distribution(command_line_options.key_distribution)),
          max_key(command_line_options.max_key),
          seed(command_line_options.seed),
          result_file(command_line_options.result_file),
          pq_settings(command_line_options.pq_settings)
#ifdef QUALITY
          ,
          histogram_file(command_line_options.histogram_file),
          log_file(command_line_options.log_file)
#endif
#ifdef WITH_PAPI
          ,
          enable_performance_counter(command_line_options.enable_performance_counter)
#endif
    {
        if (num_threads <= 0) {
            std::cerr << "Error: Number of threads must be greater than 0\n";
            throw std::invalid_argument("Invalid number of threads");
        }
        if (max_key <= 0) {
            throw std::invalid_argument("Max key must be greater than 0");
        }
        if (work_mode == WorkMode::Split) {
            if ((num_push_threads < 0 || num_push_threads > num_threads) ||
                (num_push_threads == 0 && elements_per_thread > 0)) {
                throw std::invalid_argument("Invalid number of push threads");
            }
        }
#ifdef QUALITY
        if (histogram_file.empty()) {
            throw std::invalid_argument("Histogram file must not be empty");
        }
#endif
    }

    static WorkMode parse_work_mode(char c) {
        switch (c) {
            case 'm':
                return WorkMode::Mixed;
            case 's':
                return WorkMode::Split;
            case 'i':
                return WorkMode::Increment;
            default:
                throw std::invalid_argument("Invalid work mode");
        }
    }

    static KeyDistribution parse_key_distribution(char c) {
        switch (c) {
            case 'u':
                return KeyDistribution::Uniform;
            case 'a':
                return KeyDistribution::Ascending;
            case 'd':
                return KeyDistribution::Descending;
            default:
                throw std::invalid_argument("Invalid key distribution");
        }
    }
};

void print_settings(Settings const& settings) {
    std::cout << "Threads: " << settings.num_threads << '\n'
              << "Prefill per thread: " << settings.prefill_per_thread << '\n'
              << "Elements per thread: " << settings.elements_per_thread << '\n'
              << "Operation mode: " << work_mode_name(settings.work_mode) << '\n';
    if (settings.work_mode == WorkMode::Split) {
        std::cout << "Pushing threads: " << settings.num_push_threads << '\n';
    }
    std::cout << "Key distribution: " << key_distribution_name(settings.key_distribution) << '\n'
              << "Max key: " << static_cast<unsigned long>(settings.max_key) << '\n'
              << "Seed: " << settings.seed << '\n';
#ifdef WITH_PAPI
    std::cout << "Perf. counter: " << (settings.enable_performance_counter ? "enabled" : "disabled") << '\n';
#else
    std::cout << "Perf. counter: not supported\n";
#endif
    if (!settings.result_file.empty()) {
        std::cout << "Result file: " << settings.result_file << '\n';
    }
#ifdef QUALITY
    std::cout << "Histogram file: " << settings.histogram_file << '\n';
    if (!settings.log_file.empty()) {
        std::cout << "Log file: " << settings.log_file << '\n';
    }
#else
#ifdef USE_UINT8
    std::cout << "Data type: uint8_t\n";
#else
    std::cout << "Data type: unsigned long\n";
#endif
#endif
}

#ifdef QUALITY
std::uint64_t get_tick() noexcept {
    timespec ts;
    clock_gettime(CLOCK_MONOTONIC_RAW, &ts);
    return static_cast<std::uint64_t>(ts.tv_sec) * 1000000000 + static_cast<std::uint64_t>(ts.tv_nsec);
}
#endif
#ifdef WITH_PAPI
// These are the event names used for the papi library to measure cache misses
// As these are hardware specific, you need to modify them or use generic PAPI events
static constexpr auto papi_events = std::array{"perf_raw::rc860", "perf_raw::r0864"};

void init_papi() {
    if (int ret = PAPI_library_init(PAPI_VER_CURRENT); ret != PAPI_VER_CURRENT) {
        throw std::runtime_error("Failed to initialize library");
    }
    if (int ret = PAPI_thread_init((unsigned long (*)())(pthread_self)); ret != PAPI_OK) {
        throw std::runtime_error("Failed to initialize thread support");
    }
    for (auto const& event_name : papi_events) {
        if (int ret = PAPI_query_named_event(event_name); ret != PAPI_OK) {
            throw std::runtime_error("Failed to get event code for event '" + std::string(event_name) + '\'');
        }
    }
}

void start_papi(int& event_set) {
    if (int ret = PAPI_register_thread(); ret != PAPI_OK) {
        throw std::runtime_error("Failed to register thread for PAPI");
    }
    if (int ret = PAPI_create_eventset(&event_set); ret != PAPI_OK) {
        throw std::runtime_error("Failed to create PAPI event set");
    }
    int event_code{};
    for (auto const& event_name : papi_events) {
        if (int ret = PAPI_event_name_to_code(event_name, &event_code); ret != PAPI_OK) {
            throw std::runtime_error("Failed to get PAPI event code for '" + std::string(event_name) + '\'');
        }
        if (int ret = PAPI_add_event(event_set, event_code); ret != PAPI_OK) {
            throw std::runtime_error("Failed to add PAPI event '" + std::string(event_name) + '\'');
        }
    }
    if (int ret = PAPI_start(event_set); ret != PAPI_OK) {
        throw std::runtime_error("Failed to start PAPI");
    }
}
#endif

struct Result {
    thread_coordination::time_result_type work_time{};
    long long num_pops{0};
    long long num_pushes{0};
    long long num_failed_pops{0};
#ifdef QUALITY
    operation_log::OperationLog log;
#endif
#ifdef WITH_PAPI
    std::array<long long, papi_events.size()> papi_counter{};
#endif
#ifdef MQ_COUNT_STATS
    multiqueue::Counters stat_counters{};
#endif
};

void write_results(std::vector<Result> const& results, std::ostream& out) {
    for (const auto& result : results) {
        // clang-format off
        out << result.work_time.start.time_since_epoch().count() << ','
            << result.work_time.end.time_since_epoch().count() << ','
            << result.num_pops << ','
            << result.num_pushes << ','
            << result.num_failed_pops << ','
#ifdef WITH_PAPI
            << result.papi_counter[0] << ','
            << result.papi_counter[1] << ','
#else
            << "0,0,"
#endif
#ifdef MQ_COUNT_STATS
            << result.stat_counters.locked_push_pq << ','
            << result.stat_counters.empty_pop_pqs << ','
            << result.stat_counters.locked_pop_pq << ','
            << result.stat_counters.stale_pop_pq << ','
#else
            << "0,0,0,0,"
#endif
            << '\n';
      //clang-format on
  }
}

struct ThreadData {
    thread_coordination::Context ctx;
    handle_type handle;
    std::vector<key_type> prefill;
    Result result;
};

void generate_keys(Settings const& settings, ThreadData& thread_data, std::vector<key_type>& keys) {
    auto id = thread_data.ctx.get_id();
    std::seed_seq seed{settings.seed, id};
    std::default_random_engine rng(seed);

    auto first_key = keys.begin() + id * static_cast<std::ptrdiff_t>(settings.elements_per_thread);
    auto last_key = first_key + static_cast<std::ptrdiff_t>(settings.elements_per_thread);
    auto key_dist = std::uniform_int_distribution<key_type>(0, settings.max_key);
    thread_data.ctx.execute_synchronized(
        [&key_dist, &rng](auto first, auto last) {
            std::generate(first, last, [&key_dist, &rng]() { return key_dist(rng); });
        },
        first_key, last_key);
    if (id == 0) {
        if (settings.key_distribution == KeyDistribution::Ascending) {
            std::sort(keys.begin(), keys.end());
        } else if (settings.key_distribution == KeyDistribution::Descending) {
            std::sort(keys.begin(), keys.end(), std::greater<>());
        }
    }

    thread_data.prefill.resize(settings.prefill_per_thread);
    std::generate(std::begin(thread_data.prefill), std::end(thread_data.prefill),
                  [&key_dist, &rng]() { return key_dist(rng); });
}

void push(ThreadData& thread_data, key_type key) {
#ifdef QUALITY
    auto value = operation_log::pack(thread_data.ctx.get_id(), thread_data.result.log.pushes.size());
#else
    auto value = key;
#endif
    thread_data.handle.push({key, value});
#ifdef QUALITY
    auto tick = get_tick();
    thread_data.result.log.pushes.push_back({tick, key});
#endif
}

bool try_pop(ThreadData& thread_data, pq_type::value_type& retval) {
#ifdef QUALITY
    auto tick = get_tick();
#endif
    if (!thread_data.handle.try_pop(retval)) {
        ++thread_data.result.num_failed_pops;
        return false;
    }
#ifdef QUALITY
    thread_data.result.log.pops.push_back({tick, retval.second});
#endif
    return true;
}

void prefill(ThreadData& thread_data) {
    if (thread_data.prefill.empty()) {
        return;
    }

#ifdef QUALITY
    thread_data.result.log.pushes.reserve(thread_data.prefill.size());
#endif

    thread_data.ctx.execute_synchronized([&thread_data]() {
        for (auto key : thread_data.prefill) {
            push(thread_data, key);
        }
    });
}

void execute_mixed(ThreadData& thread_data, std::vector<key_type> const& keys) {
#ifdef QUALITY
    thread_data.result.log.pushes.reserve(thread_data.result.log.pushes.size() + keys.size());
    thread_data.result.log.pops.reserve(keys.size());
#endif
    thread_data.result.work_time = thread_data.ctx.execute_synchronized_blockwise(keys.size(), [&](auto start, auto count) {
        for (auto i = start; i < start + count; ++i) {
            push(thread_data, keys[i]);
            pq_type::value_type retval;
            while (!try_pop(thread_data, retval)) {
            }
        }
        thread_data.result.num_pops += count;
        thread_data.result.num_pushes += count;
    });
}

void execute_split_push(ThreadData& thread_data, std::vector<key_type> const& keys) {
#ifdef QUALITY
    thread_data.result.log.pushes.reserve(thread_data.result.log.pushes.size() + keys.size());
#endif
    thread_data.result.work_time = thread_data.ctx.execute_synchronized_blockwise(keys.size(), [&](auto start, auto count) {
        for (auto i = start; i < start + count; ++i) {
            push(thread_data, keys[i]);
        }
        thread_data.result.num_pushes += count;
    });
}

void execute_split_pop(ThreadData& thread_data, std::atomic_llong& pop_count) {
#ifdef QUALITY
    thread_data.result.log.pops.reserve(pop_count.load(std::memory_order_relaxed));
#endif
    thread_data.result.work_time = thread_data.ctx.execute_synchronized([&thread_data, &pop_count]() {
        do {
            long long num_pops = 0;
            pq_type::value_type retval;
            while (try_pop(thread_data, retval)) {
                ++num_pops;
            }
            ++thread_data.result.num_failed_pops;
            if (num_pops == 0) {
                if (pop_count.load(std::memory_order_relaxed) <= 0) {
                    break;
                }
            } else {
                thread_data.result.num_pops += num_pops;
                if (pop_count.fetch_sub(num_pops, std::memory_order_relaxed) - num_pops <= 0) {
                    break;
                }
            }
        } while (true);
    });
}

void execute_increment(ThreadData& thread_data, std::vector<key_type> const& keys) {
#ifdef QUALITY
    thread_data.result.log.pushes.reserve(thread_data.result.log.pushes.size() + keys.size());
    thread_data.result.log.pops.reserve(keys.size());
#endif
    thread_data.result.work_time =
        thread_data.ctx.execute_synchronized_blockwise(keys.size(), [&](auto start_index, auto count) {
            for (auto i = start_index; i < start_index + count; ++i) {
                pq_type::value_type retval;
                while (!try_pop(thread_data, retval)) {
                }
                push(thread_data, retval.first + keys[i]);
            }
            thread_data.result.num_pops += count;
            thread_data.result.num_pushes += count;
        });
}

void execute_benchmark(Settings const& settings, ThreadData& thread_data, std::vector<key_type> const& keys) {
#ifdef MQ_COUNT_STATS
    thread_data.handle.reset_counters();
#endif
#ifdef WITH_PAPI
    int event_set = PAPI_NULL;
    bool papi_started = false;
    if (settings.enable_performance_counter) {
        try {
            start_papi(event_set);
            papi_started = true;
        } catch (std::exception const& e) {
            thread_data.ctx.write(std::cerr) << "Error: " << e.what() << '\n';
            thread_data.result.papi_counter.fill(-1);
        }
    }
#endif

    switch (settings.work_mode) {
        case WorkMode::Mixed:
            execute_mixed(thread_data, keys);
            break;
        case WorkMode::Split:
            if (thread_data.ctx.get_id() < settings.num_push_threads) {
                execute_split_push(thread_data, keys);
            } else {
                static std::atomic_llong pop_count =
                    static_cast<long long>(settings.prefill_per_thread + settings.elements_per_thread) *
                    settings.num_threads;
                execute_split_pop(thread_data, pop_count);
            }
            break;
        case WorkMode::Increment:
            execute_increment(thread_data, keys);
            break;
    }

#ifdef WITH_PAPI
    if (settings.enable_performance_counter && papi_started) {
        if (int ret = PAPI_stop(event_set, thread_data.result.papi_counter.data()); ret != PAPI_OK) {
            thread_data.ctx.write(std::cerr) << "Error: Failed to stop performance counters\n";
            thread_data.result.papi_counter.fill(-1);
        }
    }
#endif
#ifdef MQ_COUNT_STATS
    thread_data.result.stat_counters = thread_data.handle.get_counters();
#endif
}

void benchmark_thread(Settings const& settings, ThreadData& thread_data, std::vector<key_type>& keys) {
    if (thread_data.ctx.get_id() == 0) {
        std::clog << "Generating keys..." << std::flush;
        auto t_start = std::chrono::steady_clock::now();
        generate_keys(settings, thread_data, keys);
        auto t_end = std::chrono::steady_clock::now();
        std::clog << "done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count() << "ms)"
                  << std::endl;
        std::clog << "Prefilling... " << std::flush;
        t_start = std::chrono::steady_clock::now();
        prefill(thread_data);
        t_end = std::chrono::steady_clock::now();
        std::clog << "done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count() << "ms)"
                  << std::endl;
        std::clog << "Running benchmark..." << std::flush;
        t_start = std::chrono::steady_clock::now();
        execute_benchmark(settings, thread_data, keys);
        t_end = std::chrono::steady_clock::now();
        std::clog << "done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count() << "ms)"
                  << std::endl;
    } else {
        generate_keys(settings, thread_data, keys);
        prefill(thread_data);
        execute_benchmark(settings, thread_data, keys);
    }
}

std::vector<Result> run_benchmark(Settings const& settings) {
    auto pq = pq_type(settings.num_threads,
                      settings.prefill_per_thread * static_cast<std::size_t>(settings.num_threads), settings.pq_settings);
    std::cout << "Priority queue: ";
    pq.describe(std::cout) << '\n';
    auto keys = std::vector<key_type>(static_cast<std::size_t>(settings.num_threads) * settings.elements_per_thread);
    std::vector<Result> results(static_cast<std::size_t>(settings.num_threads));
    thread_coordination::TaskHandle task_handle(settings.num_threads, [&](auto ctx) {
        ThreadData thread_data{ctx, pq.get_handle(), {}, {}};
        benchmark_thread(settings, thread_data, keys);
        results[ctx.get_id()] = std::move(thread_data.result);
    });
    task_handle.wait();
    return results;
}

int main(int argc, char* argv[]) {
    print_header();
    CommandLineOptions options;
    cxxopts::Options cmd(argv[0]);
    cmd.add_options()
        // clang-format off
        ("h,help", "Print this help", cxxopts::value<bool>(options.help))
        ("j,threads", "The number of threads", cxxopts::value<int>(options.num_threads), "NUMBER")
        ("p,prefill", "The prefill per thread", cxxopts::value<std::size_t>(options.prefill_per_thread), "NUMBER")
        ("n,keys", "The number of keys per thread", cxxopts::value<std::size_t>(options.elements_per_thread), "NUMBER")
        ("w,work-mode", "Specify the work mode ([m]ixed, [s]plit, [i]ncrement)", cxxopts::value<char>(options.work_mode), "STRING")
        ("i,push-threads", "The number of pushing threads in split mode", cxxopts::value<int>(options.num_push_threads), "NUMBER")
        ("e,element-distribution", "Specify the element distribution ([u]niform, [a]scending, [d]escending)", cxxopts::value<char>(options.key_distribution), "STRING")
        ("m,max", "Specify the max key", cxxopts::value<key_type>(options.max_key), "NUMBER")
        ("s,seed", "Specify the initial seed", cxxopts::value<int>(options.seed), "NUMBER")
        ("x,result-file", "Path to write per-thread results to in CSV format", cxxopts::value<std::filesystem::path>(options.result_file), "PATH")
#ifdef QUALITY
        ("histogram-file", "Path to write the histogram to", cxxopts::value<std::filesystem::path>(options.histogram_file), "PATH")
        ("l,log-file", "Path to write the operation log to", cxxopts::value<std::filesystem::path>(options.log_file), "PATH")
#endif
#ifdef WITH_PAPI
        ("r,pc", "Enable performance counter", cxxopts::value<>(options.enable_performance_counter))
#endif
        // clang-format on
        ;
    pq_type::add_options(cmd, options.pq_settings);
#ifdef QUALITY
    cmd.parse_positional({"histogram-file"});
#endif
    std::optional<Settings const> settings;
    try {
        auto args = cmd.parse(argc, argv);
    } catch (std::exception const& e) {
        std::cerr << "Error parsing command line: " << e.what() << std::endl;
        std::cerr << cmd.help() << std::endl;
        return EXIT_FAILURE;
    }
    if (options.help) {
        std::clog << cmd.help() << std::endl;
        return EXIT_SUCCESS;
    }
    try {
        settings.emplace(options);
    } catch (std::exception const& e) {
        std::cerr << "Invalid options: " << e.what() << std::endl;
        std::cerr << cmd.help() << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << '\n';
    std::cout << "Command line:";
    for (int i = 0; i < argc; ++i) {
        std::cout << ' ' << argv[i];
    }
    std::cout << '\n';
    print_settings(*settings);

#ifdef WITH_PAPI
    if (settings->enable_performance_counter) {
        try {
            init_papi();
        } catch (std::exception const& e) {
            std::cerr << "Error initializing PAPI: " << e.what() << std::endl;
            return EXIT_FAILURE;
        }
    }
#endif

    std::ofstream result_file;
    if (!settings->result_file.empty()) {
        result_file = std::ofstream(settings->result_file);
        if (!result_file) {
            std::cerr << "Error: Could not open result file " << settings->result_file << " for writing" << std::endl;
            return EXIT_FAILURE;
        }
    }

#ifdef QUALITY
    std::ofstream log_file;
    if (!settings->log_file.empty()) {
        log_file = std::ofstream(settings->log_file);
        if (!log_file) {
            std::cerr << "Error: Could not open log file " << settings->log_file << " for writing" << std::endl;
            return EXIT_FAILURE;
        }
    }
    auto histogram_file = std::ofstream(settings->histogram_file);
    if (!histogram_file) {
        std::cerr << "Error: Could not open histogram file " << settings->histogram_file << " for writing" << std::endl;
        return EXIT_FAILURE;
    }

#endif

    auto results = run_benchmark(*settings);

    if (!settings->result_file.empty()) {
        write_results(results, result_file);
    }

#ifdef QUALITY
    std::vector<operation_log::OperationLog> logs;
    for (auto& result : results) {
        logs.push_back(std::move(result.log));
    }
    if (!settings->log_file.empty()) {
        std::clog << "Writing operation logs..." << std::flush;
        auto t_start = std::chrono::steady_clock::now();
        operation_log::write_logs(logs, log_file);
        log_file.close();
        auto t_end = std::chrono::steady_clock::now();
        std::clog << "done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count()
                  << "ms)\n";
    }
    try {
        operation_log::verify_logs(logs);
    } catch (std::exception const& e) {
        std::cerr << "Error verifying logs: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
    std::clog << "Computing quality histogram..." << std::flush;
    auto t_start = std::chrono::steady_clock::now();
    auto histogram = operation_log::to_histogram(logs);
    auto t_end = std::chrono::steady_clock::now();
    std::clog << "done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count() << "ms)\n";
    operation_log::write_histogram(histogram, histogram_file);
#endif
}
