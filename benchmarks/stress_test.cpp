#include "build_info.hpp"
#include "priority_queue_selection.hpp"
#include "thread_coordination.hpp"

#include "cxxopts.hpp"

#ifdef QUALITY
#undef USE_UINT8
#include "operation_log.hpp"
#else
#include "counting_handle.hpp"
#ifdef WITH_PAPI
#define USE_PAPI
#include <papi.h>
#include <pthread.h>
#endif
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

#ifdef USE_UINT8
using key_type = std::uint8_t;
using value_type = std::uint8_t;
#else
using key_type = unsigned long;
using value_type = unsigned long;
#endif

using pq_type = PriorityQueue<key_type, value_type, true>;
#ifdef QUALITY
using handle_type = operation_log::LoggingHandle<pq_type>;
#else
using handle_type = CountingHandle<pq_type>;
#endif
using ctx_type = thread_coordination::Context;

enum class WorkMode { Mixed, Push, Pop, Increment };
enum class KeyDistribution { Uniform, Ascending, Descending };

[[nodiscard]] auto work_mode_name(WorkMode work_mode) {
    switch (work_mode) {
        case WorkMode::Mixed:
            return "mixed";
        case WorkMode::Push:
            return "push";
        case WorkMode::Pop:
            return "pop";
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
    int num_threads = 4;
    std::size_t prefill_per_thread = 1 << 20;
    std::size_t elements_per_thread = 1 << 24;
    WorkMode work_mode = WorkMode::Mixed;
    KeyDistribution key_distribution = KeyDistribution::Uniform;
    pq_type::config_type pq_settings{};
#ifdef USE_UINT8
    key_type max_key = std::numeric_limits<key_type>::max();
#else
    key_type max_key = 1 << 30;
#endif
    int seed = 1;
#ifdef QUALITY
    std::filesystem::path log_file;
#endif
#ifdef USE_PAPI
    std::vector<std::string> papi_events;
#endif
};

bool validate_settings(Settings const& settings) {
    if (settings.num_threads <= 0) {
        std::cerr << "Error: Number of threads must be greater than 0\n";
        return false;
    }
    if (settings.max_key <= 0) {
        std::cerr << "Error: Max key must be greater than 0\n";
        return false;
    }
#ifdef USE_PAPI
    for (auto const& name : settings.papi_events) {
        if (PAPI_query_named_event(name.c_str()) != PAPI_OK) {
            std::cerr << "Error: Invalid PAPI event '" << name << "'\n";
            return false;
        }
    }
#endif
    return true;
}

void write_settings(Settings const& settings, std::ostream& out) {
    out << "Threads: " << settings.num_threads << '\n'
        << "Prefill per thread: " << settings.prefill_per_thread << '\n'
        << "Elements per thread: " << settings.elements_per_thread << '\n'
        << "Operation mode: " << work_mode_name(settings.work_mode) << '\n';
    out << "Key distribution: " << key_distribution_name(settings.key_distribution) << '\n'
        << "Max key: " << static_cast<unsigned long>(settings.max_key) << '\n'
        << "Seed: " << settings.seed << '\n';
#ifdef USE_PAPI
    if (!settings.papi_events.empty()) {
        out << "PAPI events: ";
        std::copy(settings.papi_events.begin(), settings.papi_events.end(),
                  std::ostream_iterator<std::string>(out, " "));
        out << '\n';
    }
#endif
#ifdef QUALITY
    if (!settings.log_file.empty()) {
        out << "Log operations to: " << settings.log_file << '\n';
    }
#else
#ifdef USE_UINT8
    out << "Data type: uint8_t\n";
#else
    out << "Data type: unsigned long\n";
#endif
#endif
}

#ifdef USE_PAPI
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

void prepare(Settings const& settings, ctx_type const& ctx, handle_type& handle, std::vector<key_type>& keys) {
    auto id = ctx.get_id();
    std::seed_seq seed{settings.seed, id};
    std::default_random_engine rng(seed);

    auto first_key = keys.begin() + id * static_cast<std::ptrdiff_t>(settings.elements_per_thread);
    auto last_key = first_key + static_cast<std::ptrdiff_t>(settings.elements_per_thread);
    auto key_dist = std::uniform_int_distribution<key_type>(0, settings.max_key);
    ctx.execute_synchronized(
        [&](auto first, auto last) {
            switch (settings.key_distribution) {
                case KeyDistribution::Uniform: {
                    std::generate(first, last, [&key_dist, &rng]() { return key_dist(rng); });
                    break;
                }
                case KeyDistribution::Ascending: {
                    auto i = static_cast<std::size_t>(std::distance(keys.begin(), first));
                    for (; first != last; ++first, ++i) {
                        *first = static_cast<key_type>((i * settings.max_key) / keys.size() + 1);
                    }
                    break;
                }
                case KeyDistribution::Descending: {
                    auto i = static_cast<std::size_t>(std::distance(first, keys.end()) - 1);
                    for (; first != last; ++first, --i) {
                        *first = static_cast<key_type>((i * settings.max_key) / keys.size() + 1);
                    }
                    break;
                }
            }
        },
        first_key, last_key);

    std::vector<key_type> prefill_keys(settings.prefill_per_thread);
    std::generate(std::begin(prefill_keys), std::end(prefill_keys), [&key_dist, &rng]() { return key_dist(rng); });
    ctx.execute_synchronized([&]() {
        for (auto const& key : prefill_keys) {
            handle.push(key);
        }
    });
}

auto execute_mixed(ctx_type const& ctx, handle_type& handle, std::vector<key_type> const& keys) {
    return ctx.execute_synchronized_blockwise(keys.begin(), keys.end(), [&](auto first, auto last) {
        for (; first != last; ++first) {
            handle.push(*first);
            while (!handle.try_pop()) {
            }
        }
    });
}

auto execute_push(ctx_type const& ctx, handle_type& handle, std::vector<key_type> const& keys) {
    return ctx.execute_synchronized_blockwise(keys.begin(), keys.end(), [&](auto first, auto last) {
        for (; first != last; ++first) {
            handle.push(*first);
        }
    });
}

auto execute_pop(ctx_type const& ctx, handle_type& handle, std::vector<key_type> const& keys) {
    return ctx.execute_synchronized_blockwise(keys.begin(), keys.end(), [&](auto first, auto last) {
        for (; first != last; ++first) {
            while (!handle.try_pop()) {
            }
        }
    });
}

auto execute_increment(ctx_type const& ctx, handle_type& handle, std::vector<key_type> const& keys) {
    return ctx.execute_synchronized_blockwise(keys.begin(), keys.end(), [&](auto first, auto last) {
        for (; first != last; ++first) {
            auto retval = handle.try_pop();
            while (!retval) {
                retval = handle.try_pop();
            }
            handle.push(retval->first + *first);
        }
    });
}

struct BenchmarkResult {
    thread_coordination::time_result_type work_time{};
#ifdef USE_PAPI
    std::vector<long long> papi_counter{};
#endif
#ifdef MQ_COUNT_STATS
    multiqueue::Counters mq_stats{};
#endif
#ifdef QUALITY
    operation_log::OperationLog op_log{};
#else
    OperationCounters op_counters{};
#endif
};

#ifndef QUALITY
void write_result(Settings const& settings, std::vector<BenchmarkResult> const& stats, std::ostream& out) {
    out << "thread,start_time,end_time,pushes,pops,failed_pops,locked_push_pq,empty_pop_pqs,"
           "locked_pop_pq,stale_pop_pq";
#ifdef USE_PAPI
    for (auto const& e : settings.papi_events) {
        out << ',' << e;
    }
#endif
    out << '\n';
    for (int i = 0; i < settings.num_threads; ++i) {
        auto const& s = stats[static_cast<std::size_t>(i)];
        // clang-format off
        out << i << ','
            << s.work_time.start.time_since_epoch().count() << ','
            << s.work_time.end.time_since_epoch().count() << ','
            << s.op_counters.push << ','
            << s.op_counters.pop << ','
            << s.op_counters.failed_pop;
        // clang-format on
#ifdef MQ_COUNT_STATS
        // clang-format off
        out << ','
            << s.mq_stats.locked_push_pq << ','
            << s.mq_stats.empty_pop_pqs << ','
            << s.mq_stats.locked_pop_pq << ','
            << s.mq_stats.stale_pop_pq;
        // clang-format on
#else
        out << ",0,0,0,0";
#endif
#ifdef USE_PAPI
        for (auto c : s.papi_counter) {
            out << ',' << c;
        }
#endif
        out << '\n';
    }
}
#endif

BenchmarkResult execute_benchmark(Settings const& settings, ctx_type const& ctx, handle_type& handle,
                                  std::vector<key_type> const& keys) {
    BenchmarkResult result;
#ifdef MQ_COUNT_STATS
    handle.reset_counters();
#endif
#ifdef USE_PAPI
    int event_set = PAPI_NULL;
    bool papi_started = false;
    if (!settings.papi_events.empty()) {
        result.papi_counter.resize(settings.papi_events.size());
        try {
            event_set = start_papi(settings.papi_events);
            papi_started = true;
        } catch (std::exception const& e) {
            ctx.write(std::cerr) << "Error: " << e.what() << '\n';
            std::fill(result.papi_counter.begin(), result.papi_counter.end(), -1);
        }
    }
#endif

    auto benchmark = [&settings]() {
        switch (settings.work_mode) {
            case WorkMode::Mixed:
                return &execute_mixed;
            case WorkMode::Push:
                return &execute_push;
            case WorkMode::Pop:
                return &execute_pop;
            case WorkMode::Increment:
                return &execute_increment;
        }
        __builtin_unreachable();
    }();

    result.work_time = benchmark(ctx, handle, keys);

#ifdef USE_PAPI
    if (papi_started) {
        if (int ret = PAPI_stop(event_set, result.papi_counter.data()); ret != PAPI_OK) {
            ctx.write(std::cerr) << "Error: Failed to stop performance counters\n";
            std::fill(result.papi_counter.begin(), result.papi_counter.end(), -1);
        }
    }
#endif
#ifdef MQ_COUNT_STATS
    result.mq_stats = handle.get_counters();
#endif
#ifdef QUALITY
    result.op_log = handle.extract_log();
#else
    result.op_counters = handle.get_operation_counts();
#endif
    return result;
}

std::vector<BenchmarkResult> run_benchmark(Settings const& settings) {
    auto pq =
        pq_type(settings.num_threads, settings.prefill_per_thread * static_cast<std::size_t>(settings.num_threads),
                settings.pq_settings);
    std::clog << "Priority queue: ";
    pq.describe(std::clog) << '\n';
    auto keys = std::vector<key_type>(static_cast<std::size_t>(settings.num_threads) * settings.elements_per_thread);
    std::vector<BenchmarkResult> results(static_cast<std::size_t>(settings.num_threads));
    thread_coordination::TaskHandle task_handle(settings.num_threads, [&](auto ctx) {
#ifdef QUALITY
        auto handle = handle_type(ctx.get_id(), pq);
        handle.reserve_push_log(2 * (settings.prefill_per_thread + settings.elements_per_thread));
        handle.reserve_pop_log(2 * settings.elements_per_thread);
#else
        auto handle = handle_type(pq);
#endif
        if (ctx.get_id() == 0) {
            std::clog << "Preparing..." << std::flush;
            auto t_start = std::chrono::steady_clock::now();
            prepare(settings, ctx, handle, keys);
            auto t_end = std::chrono::steady_clock::now();
            std::clog << "done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count()
                      << "ms)" << std::endl;
            std::clog << "Running benchmark..." << std::flush;
            t_start = std::chrono::steady_clock::now();
            results[static_cast<std::size_t>(ctx.get_id())] = execute_benchmark(settings, ctx, handle, keys);
            t_end = std::chrono::steady_clock::now();
            std::clog << "done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count()
                      << "ms)" << std::endl;
        } else {
            prepare(settings, ctx, handle, keys);
            results[static_cast<std::size_t>(ctx.get_id())] = execute_benchmark(settings, ctx, handle, keys);
        }
    });
    task_handle.wait();
    return results;
}

WorkMode parse_work_mode(char c) {
    switch (c) {
        case 'm':
            return WorkMode::Mixed;
        case 'u':
            return WorkMode::Push;
        case 'o':
            return WorkMode::Pop;
        case 'i':
            return WorkMode::Increment;
        default:
            throw std::invalid_argument("Invalid work mode");
    }
}

KeyDistribution parse_key_distribution(char c) {
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

int main(int argc, char* argv[]) {
    write_build_info(std::clog);
    Settings settings;
    cxxopts::Options cmd(argv[0]);
    cmd.add_options()
        // clang-format off
        ("h,help", "Print this help", cxxopts::value<bool>())
        ("j,threads", "The number of threads", cxxopts::value<int>(settings.num_threads), "NUMBER")
        ("p,prefill", "The prefill per thread", cxxopts::value<std::size_t>(settings.prefill_per_thread), "NUMBER")
        ("n,keys", "The number of keys per thread", cxxopts::value<std::size_t>(settings.elements_per_thread), "NUMBER")
        ("w,work-mode", "Specify the work mode ([m]ixed, p[u]sh, p[o]p, [i]ncrement)", cxxopts::value<char>(), "STRING")
        ("d,key-dist", "Specify the key distribution ([u]niform, [a]scending, [d]escending)", cxxopts::value<char>(), "STRING")
        ("m,max", "Specify the max key", cxxopts::value<key_type>(settings.max_key), "NUMBER")
        ("s,seed", "Specify the initial seed", cxxopts::value<int>(settings.seed), "NUMBER")
#ifdef QUALITY
        ("l,log-file", "Path to write the operation log to", cxxopts::value<std::filesystem::path>(settings.log_file), "PATH")
#endif
#ifdef USE_PAPI
        ("r,pc", "Performance counters", cxxopts::value<std::vector<std::string>>(settings.papi_events))
#endif
        // clang-format on
        ;
    pq_type::add_options(cmd, settings.pq_settings);
    try {
        auto args = cmd.parse(argc, argv);
        if (args.count("help") > 0) {
            std::clog << cmd.help() << std::endl;
            return EXIT_SUCCESS;
        }
        if (args.count("work-mode") > 0) {
            settings.work_mode = parse_work_mode(args["work-mode"].as<char>());
        }
        if (args.count("key-dist") > 0) {
            settings.key_distribution = parse_key_distribution(args["key-dist"].as<char>());
        }
    } catch (std::exception const& e) {
        std::cerr << "Error parsing command line: " << e.what() << std::endl;
        std::cerr << cmd.help() << std::endl;
        return EXIT_FAILURE;
    }

    std::clog << '\n';
    std::clog << "Command line:";
    for (int i = 0; i < argc; ++i) {
        std::clog << ' ' << argv[i];
    }
    std::clog << '\n';
    write_settings(settings, std::clog);

#ifdef USE_PAPI
    if (!settings.papi_events.empty()) {
        if (int ret = PAPI_library_init(PAPI_VER_CURRENT); ret != PAPI_VER_CURRENT) {
            std::cerr << "Error: Failed to initialize PAPI library\n";
            return EXIT_FAILURE;
        }
        if (int ret = PAPI_thread_init((unsigned long (*)())(pthread_self)); ret != PAPI_OK) {
            std::cerr << "Error: Failed to initialize PAPI thread support\n";
            return EXIT_FAILURE;
        }
    }
#endif

    if (!validate_settings(settings)) {
        return EXIT_FAILURE;
    }

#ifdef QUALITY
    std::ofstream log_file;
    if (!settings.log_file.empty()) {
        log_file = std::ofstream(settings.log_file);
        if (!log_file) {
            std::cerr << "Error: Could not open log file " << settings.log_file << " for writing" << std::endl;
            return EXIT_FAILURE;
        }
    }
#endif

    auto results = run_benchmark(settings);

#ifdef QUALITY
    std::vector<operation_log::OperationLog> logs;
    for (auto&& result : results) {
        logs.push_back(std::move(result.op_log));
    }
    if (!settings.log_file.empty()) {
        std::clog << "Writing operation logs..." << std::flush;
        auto t_start = std::chrono::steady_clock::now();
        operation_log::write_logs(logs, log_file);
        log_file.close();
        auto t_end = std::chrono::steady_clock::now();
        std::clog << "done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count()
                  << "ms)\n";
    }
    if (!operation_log::verify_logs(logs)) {
        return false;
    }
    std::clog << "Replaying logs..." << std::flush;
    auto t_start = std::chrono::steady_clock::now();
    auto metrics = operation_log::replay_logs(logs);
    auto t_end = std::chrono::steady_clock::now();
    std::clog << "done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count() << "ms)\n";
    operation_log::write_metrics(metrics, std::cout);
#else
    write_result(settings, results, std::cout);
#endif
}
