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
#ifdef QUALITY
#include <filesystem>
#endif

void print_header() {
    std::clog << "Built on " << __DATE__ << ' ' << __TIME__ << " with:\n";
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
#ifdef WITH_PAPI
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

using thread_coordination::timepoint_type;

enum class WorkMode { Mixed, Split, Increment };
enum class KeyDistribution { Uniform, Ascending, Descending };

struct Settings {
    int num_threads = 4;
    std::size_t prefill_per_thread = 1 << 20;
    std::size_t elements_per_thread = 1 << 24;
    WorkMode work_mode = WorkMode::Mixed;
    int num_push_threads = 1;
    KeyDistribution key_distribution = KeyDistribution::Uniform;
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

    bool set_work_mode(char c) {
        switch (c) {
            case 'm':
                work_mode = WorkMode::Mixed;
                return true;
            case 's':
                work_mode = WorkMode::Split;
                return true;
            case 'i':
                work_mode = WorkMode::Increment;
                return true;
        }
        return false;
    }

    bool set_key_distribution(char c) {
        switch (c) {
            case 'u':
                key_distribution = KeyDistribution::Uniform;
                return true;
            case 'a':
                key_distribution = KeyDistribution::Ascending;
                return true;
            case 'd':
                key_distribution = KeyDistribution::Descending;
                return true;
        }
        return false;
    }

    [[nodiscard]] auto work_mode_str() const {
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
    [[nodiscard]] auto key_distribution_str() const {
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

    [[nodiscard]] bool validate() const {
        if (num_threads <= 0) {
            std::cerr << "Error: Number of threads must be greater than 0\n";
            return false;
        }
        if (max_key == 0) {
            std::cerr << "Error: Max key must be greater than 0\n";
            return false;
        }
        if (work_mode == WorkMode::Split) {
            if ((num_push_threads < 0 || num_push_threads > num_threads) ||
                (num_push_threads == 0 && elements_per_thread > 0)) {
                std::cerr << "Error: Invalid number of push threads\n";
                return false;
            }
            if (num_push_threads == 0 && elements_per_thread > 0) {
                std::cerr << "Error: Invalid number of push threads\n";
                return false;
            }
        }
        return true;
    }
};

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
        throw std::runtime_error("Failed to initialize PAPI");
    }
    if (int ret = PAPI_thread_init((unsigned long (*)())(pthread_self)); ret != PAPI_OK) {
        throw std::runtime_error("Failed to initialize PAPI thread support");
    }
    for (auto const& event_name : papi_events) {
        if (int ret = PAPI_query_named_event(event_name); ret != PAPI_OK) {
            throw std::runtime_error("Failed to get PAPI event code for event '" + std::string(event_name) + '\'');
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

void print_settings(Settings const& settings) {
    std::cout << "Threads: " << settings.num_threads << '\n'
              << "Prefill per thread: " << settings.prefill_per_thread << '\n'
              << "Elements per thread: " << settings.elements_per_thread << '\n'
              << "Operation mode: " << settings.work_mode_str();
    if (settings.work_mode == WorkMode::Split) {
        std::cout << " (" << settings.num_push_threads << " push)";
    }
    std::cout << '\n';
    std::cout << "Key distribution: " << settings.key_distribution_str() << '\n'
              << "Max key: " << static_cast<unsigned long>(settings.max_key) << '\n'
              << "Seed: " << settings.seed << '\n';
#ifdef WITH_PAPI
    std::cout << "Perf. counter: " << (settings.enable_performance_counter ? "enabled" : "disabled") << '\n';
#else
    std::cout << "Perf. counter: unavailable\n";
#endif
#ifdef USE_UINT8
    std::cout << "Data type: uint8_t\n";
#else
    std::cout << "Data type: unsigned long\n";
#endif
    std::cout << '\n' << '\n';
}

struct Result {
    timepoint_type start_time{timepoint_type::max()};
    timepoint_type end_time{timepoint_type::min()};
    long long num_failed_pops{0};
    long long num_pops{0};
    long long num_pushes{0};
#ifdef WITH_PAPI
    std::array<long long, papi_events.size()> papi_counter{};
#endif
#ifdef MQ_COUNT_STATS
    multiqueue::Counters stat_counters;
#endif
#ifdef QUALITY
    operation_log::OperationLog log;
#endif
};

void print_result(std::vector<Result> const& results) {
    auto start_time = std::min_element(results.begin(), results.end(), [](auto const& lhs, auto const& rhs) {
                          return lhs.start_time < rhs.start_time;
                      })->start_time;
    auto end_time = std::max_element(results.begin(), results.end(), [](auto const& lhs, auto const& rhs) {
                        return lhs.end_time < rhs.end_time;
                    })->end_time;
    std::cout << "Work time (s): " << std::setprecision(3)
              << std::chrono::duration<double>(end_time - start_time).count() << '\n';
    /* #ifdef WITH_PAPI */
    /*     std::cout << "L1d cache misses: " << result.l1d_cache_misses << '\n'; */
    /*     std::cout << "L2 cache misses: " << result.l2_cache_misses << '\n'; */
    /* #endif */
    /*     std::cout << result.success ? "Benchmark successful" : "Benchmark failed" << '\n'; */
}

/* void write_csv(Settings const& settings, Result const& result) { */
/*     std::ofstream out(settings.result_file); */
/*     if (!out) { */
/*         throw std::runtime_error("Failed to open file: " + file.string()); */
/*     } */
/*     out << "threads,prefill,elements,work-mode,push-threads,key-distribution,max-key,seed,work-time," */
/*            "failed-pops,l1d-cache-misses,l2-cache-misses,locked-push-pq,empty_pop_pqs,locked_pop_pq,stale_pop_" */
/*            "pq"; */
/*     out << settings.num_threads << ',' << settings.prefill_per_thread << ',' << settings.elements_per_thread << ','
 */
/*         << settings.work_mode_str() << ',' << settings.num_push_threads << ',' << settings.key_distribution_str() <<
 * ',' */
/*         << static_cast<unsigned long>(settings.max_key) << ',' << settings.seed << ',' << std::fixed */
/*         << std::setprecision(3) */
/*         << std::chrono::duration<double>(result.work_time.second - result.work_time.first).count() << ',' */
/*         << result.num_failed_pops; */
/* #ifdef WITH_PAPI */
/*     if (settings.enable_performance_counter) { */
/*         out << ',' << result.l1d_cache_misses << ',' << result.l2_cache_misses; */
/*     } else { */
/*         out << ",-1,-1"; */
/*     } */
/* #else */
/*     out << ",-1,-1"; */
/* #endif */
/* #ifdef MQ_COUNT_STATS */
/*     out << ',' << result.counters.locked_push_pq << ',' << result.counters.empty_pop_pqs << ',' */
/*         << result.counters.locked_pop_pq << ',' << result.counters.stale_pop_pq; */
/* #else */
/*     out << ",-1,-1,-1,-1"; */
/* #endif */
/*     out << ',' << result.success << std::endl; */
/* } */

struct ThreadData {
    thread_coordination::Context ctx;
    handle_type handle;
    std::vector<key_type> prefill;
    Result result{};
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
    thread_data.result.log.pops.push_back({tick, retval.first, retval.second});
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
    auto [t_start, t_end] = thread_data.ctx.execute_synchronized_blockwise(keys.size(), [&](auto start, auto count) {
        for (auto i = start; i < start + count; ++i) {
            push(thread_data, keys[i]);
            pq_type::value_type retval;
            while (!try_pop(thread_data, retval)) {
            }
        }
        thread_data.result.num_pops += count;
        thread_data.result.num_pushes += count;
    });
    thread_data.result.start_time = t_start;
    thread_data.result.end_time = t_end;
}

void execute_split_push(ThreadData& thread_data, std::vector<key_type> const& keys) {
#ifdef QUALITY
    thread_data.result.log.pushes.reserve(thread_data.result.log.pushes.size() + keys.size());
#endif
    auto [t_start, t_end] = thread_data.ctx.execute_synchronized_blockwise(keys.size(), [&](auto start, auto count) {
        for (auto i = start; i < start + count; ++i) {
            push(thread_data, keys[i]);
        }
        thread_data.result.num_pushes += count;
    });
    thread_data.result.start_time = t_start;
    thread_data.result.end_time = t_end;
}

void execute_split_pop(ThreadData& thread_data, std::atomic_llong& pop_count) {
#ifdef QUALITY
    thread_data.result.log.pops.reserve(pop_count.load(std::memory_order_relaxed));
#endif
    auto [t_start, t_end] = thread_data.ctx.execute_synchronized([&thread_data, &pop_count]() {
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
    thread_data.result.start_time = t_start;
    thread_data.result.end_time = t_end;
}

void execute_increment(ThreadData& thread_data, std::vector<key_type> const& keys) {
#ifdef QUALITY
    thread_data.result.log.pushes.reserve(thread_data.result.log.pushes.size() + keys.size());
    thread_data.result.log.pops.reserve(keys.size());
#endif
    auto [t_start, t_end] =
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
    thread_data.result.start_time = t_start;
    thread_data.result.end_time = t_end;
}

void execute_benchmark(Settings const& settings, ThreadData& thread_data, std::vector<key_type> const& keys) {
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
    thread_data.result.counters = thread_data.handle.get_counters();
#endif
}

void benchmark_thread(Settings const& settings, ThreadData& thread_data, std::vector<key_type>& keys) {
    if (thread_data.ctx.get_id() == 0) {
        std::clog << "Generating keys..." << std::flush;
    }
    generate_keys(settings, thread_data, keys);
    if (thread_data.ctx.get_id() == 0) {
        std::clog << "done\nPrefilling..." << std::flush;
    }
    prefill(thread_data);
#ifdef MQ_COUNT_STATS
    thread_data.handle.reset_counters();
#endif
    if (thread_data.ctx.get_id() == 0) {
        std::clog << "done\nRunning benchmark..." << std::flush;
    }
    execute_benchmark(settings, thread_data, keys);
    if (thread_data.ctx.get_id() == 0) {
        std::clog << "done" << std::endl;
    }
}

void run_benchmark(Settings const& settings, pq_type& pq) {
    auto keys = std::vector<key_type>(static_cast<std::size_t>(settings.num_threads) * settings.elements_per_thread);
    std::vector<Result> results(settings.num_threads);
    thread_coordination::TaskHandle task_handle(settings.num_threads, [&](auto ctx) {
        ThreadData thread_data{ctx, pq.get_handle()};
        benchmark_thread(settings, thread_data, keys);
        results[ctx.get_id()] = std::move(thread_data.result);
    });
    task_handle.wait();

#ifdef QUALITY
    std::vector<operation_log::OperationLog> logs;
    for (auto const& result : results) {
        logs.push_back(std::move(result.log));
    }
    auto sorted_logs = operation_log::sort_logs(logs);
    if (!settings.log_file.empty()) {
        std::clog << "Writing operation logs..." << std::flush;
        try {
            operation_log::write_logs(sorted_logs, settings.log_file);
            std::clog << "done" << std::endl;
        } catch (std::exception const& e) {
            std::clog << "failed\n" << std::endl;
            std::cerr << "Error: " << e.what() << std::endl;
        }
    }
    std::clog << "Computing quality histogram..." << std::flush;
    try {
        auto histogram = operation_log::to_histogram(sorted_logs);
        std::clog << "done\nWriting quality histogram..." << std::flush;
        std::ofstream out(settings.histogram_file);
        if (!out) {
            throw std::runtime_error("Failed to open file: " + settings.histogram_file.string());
        }
        out << sorted_logs.pops.size() << '\n';
        for (std::size_t i = 0; i < std::max(histogram.ranks.size(), histogram.delays.size()); ++i) {
            out << i << ' ' << (i < histogram.ranks.size() ? histogram.ranks[i] : 0) << ' '
                << (i < histogram.delays.size() ? histogram.delays[i] : 0) << '\n';
        }
        std::clog << "done" << std::endl;
    } catch (std::runtime_error const& e) {
        std::clog << "failed" << std::endl;
        std::cerr << "Error: " << e.what() << std::endl;
    }
#endif
    print_result(results);

    if (!settings.result_file.empty()) {
        try {
            /* write_csv(results, settings.result_file); */
        } catch (std::exception const& e) {
            std::cerr << "Error: " << e.what() << std::endl;
        }
    }
}

int main(int argc, char* argv[]) {
    print_header();

    std::clog << "\nCommand line:";
    for (int i = 0; i < argc; ++i) {
        std::clog << ' ' << argv[i];
    }
    std::clog << '\n' << '\n';

    Settings settings;
    cxxopts::Options options(argv[0]);
    options.add_options()
        // clang-format off
        ("h,help", "Print this help")
        ("j,threads", "The number of threads", cxxopts::value<int>(settings.num_threads), "NUMBER")
        ("p,prefill", "The prefill per thread", cxxopts::value<std::size_t>(settings.prefill_per_thread), "NUMBER")
        ("n,keys", "The number of keys per thread", cxxopts::value<std::size_t>(settings.elements_per_thread), "NUMBER")
        ("w,work-mode", "Specify the work mode ([m]ixed, [s]plit, [i]ncrement)", cxxopts::value<char>(), "STRING")
        ("i,push-threads", "The number of pushing threads in split mode", cxxopts::value<int>(settings.num_push_threads), "NUMBER")
        ("e,element-distribution", "Specify the element distribution ([u]niform, [a]scending, [d]escending)", cxxopts::value<char>(), "STRING")
        ("m,max", "Specify the max key", cxxopts::value<key_type>(settings.max_key), "NUMBER")
        ("s,seed", "Specify the initial seed", cxxopts::value<int>(settings.seed), "NUMBER")
        ("x,result-file", "Path to write results to in CSV format", cxxopts::value<std::filesystem::path>(settings.result_file), "PATH")
#ifdef QUALITY
        ("q,histogram-file", "Path to write the histogram to", cxxopts::value<std::filesystem::path>(settings.histogram_file), "PATH")
        ("l,log-file", "Path to write the operation log to", cxxopts::value<std::filesystem::path>(settings.log_file), "PATH")
#endif
#ifdef WITH_PAPI
        ("r,pc", "Enable performance counter", cxxopts::value<>(settings.enable_performance_counter))
#endif
        // clang-format on
        ;
    pq_type::add_options(options);

    cxxopts::ParseResult result;
    try {
        result = options.parse(argc, argv);
    } catch (cxxopts::OptionParseException const& e) {
        std::cerr << "Error parsing arguments: " << e.what() << '\n';
        std::cerr << options.help() << std::endl;
        return 1;
    }
    if (result.count("help") > 0) {
        std::cerr << options.help() << std::endl;
        return 0;
    }
    if (result.count("work-mode") > 0) {
        auto work_mode = result["work-mode"].as<char>();
        if (!settings.set_work_mode(work_mode)) {
            std::cerr << "Invalid work mode: " << result["work-mode"].as<char>() << '\n';
            std::cerr << options.help() << std::endl;
            return 1;
        }
    }
    if (result.count("element-distribution") > 0) {
        auto key_distribution = result["element-distribution"].as<char>();
        if (!settings.set_key_distribution(key_distribution)) {
            std::cerr << "Invalid element distribution: " << result["element-distribution"].as<char>() << '\n';
            std::cerr << options.help() << std::endl;
            return 1;
        }
    }

    print_settings(settings);

    if (!settings.validate()) {
        return 1;
    }

#ifdef WITH_PAPI
    if (settings.enable_performance_counter) {
        try {
            init_papi();
        } catch (std::exception const& e) {
            std::cerr << "Error: " << e.what() << '\n';
            std::clog << "Performance counter disabled" << std::endl;
            settings.enable_performance_counter = false;
        }
    }
#endif

    auto pq = pq_type(settings.num_threads,
                      settings.prefill_per_thread * static_cast<std::size_t>(settings.num_threads), result);
    run_benchmark(settings, pq);
}
