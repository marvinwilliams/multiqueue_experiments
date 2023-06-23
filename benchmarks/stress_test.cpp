#ifdef QUALITY
#undef WITH_PAPI
#undef USE_UINT8
#undef MQ_COUNT_STATS
#define WITH_OPERATION_LOG
#endif

#include "common.hpp"
#include "priority_queue_selection.hpp"
#include "thread_coordination.hpp"

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

#ifdef USE_UINT8
using key_type = std::uint8_t;
using value_type = std::uint8_t;
#else
using key_type = unsigned long;
using value_type = unsigned long;
#endif

using pq_type = PriorityQueue<key_type, value_type, true>;
#ifdef QUALITY
using handle_type = Handle<pq_type, true>;
#else
using handle_type = Handle<pq_type, false>;
#endif
using thread_context_type = thread_coordination::Context;

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
    int num_threads = 4;
    std::size_t prefill_per_thread = 1 << 20;
    std::size_t elements_per_thread = 1 << 24;
    WorkMode work_mode = WorkMode::Mixed;
    int num_push_threads = 1;
    KeyDistribution key_distribution = KeyDistribution::Uniform;
    pq_type::config_type pq_settings{};
#ifdef USE_UINT8
    key_type max_key = std::numeric_limits<key_type>::max();
#else
    key_type max_key = 1 << 30;
#endif
    int seed = 1;
    std::filesystem::path result_file;
#ifdef QUALITY
    std::filesystem::path log_file;
#endif
#ifdef WITH_PAPI
    std::vector<std::string> papi_events;
#endif
};

bool validate_settings(Settings const& settings) {
    if (settings.num_threads <= 0) {
        std::cerr << "Error: Number of threads must be greater than 0\n";
        return false;
    }
    if (settings.max_key <= 0) {
        std::cerr << "Max key must be greater than 0\n";
        return false;
    }
    if (settings.work_mode == WorkMode::Split) {
        if ((settings.num_push_threads < 0 || settings.num_push_threads > settings.num_threads) ||
            (settings.num_push_threads == 0 && settings.elements_per_thread > 0)) {
            std::cerr << "Invalid number of push threads";
            return false;
        }
    }
    if (settings.result_file.empty()) {
        std::cerr << "Result file must not be empty\n";
        return false;
    }
#ifdef WITH_PAPI
    for (auto const& name : settings.papi_events) {
        if (PAPI_query_named_event(name.c_str()) != PAPI_OK) {
            std::cerr << "Invalid PAPI event: " << name << '\n';
            return false;
        }
    }
#endif
    return true;
}

WorkMode parse_work_mode(char c) {
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

void write_settings(Settings const& settings, std::ostream& out) {
    out << "Threads: " << settings.num_threads << '\n'
        << "Prefill per thread: " << settings.prefill_per_thread << '\n'
        << "Elements per thread: " << settings.elements_per_thread << '\n'
        << "Operation mode: " << work_mode_name(settings.work_mode) << '\n';
    if (settings.work_mode == WorkMode::Split) {
        out << "Pushing threads: " << settings.num_push_threads << '\n';
    }
    out << "Key distribution: " << key_distribution_name(settings.key_distribution) << '\n'
        << "Max key: " << static_cast<unsigned long>(settings.max_key) << '\n'
        << "Seed: " << settings.seed << '\n';
#ifdef WITH_PAPI
    if (!settings.papi_events.empty()) {
        out << "PAPI events: ";
        std::copy(settings.papi_events.begin(), settings.papi_events.end(),
                  std::ostream_iterator<std::string>(out, " "));
        out << '\n';
    }
#endif
    out << "Result file: " << settings.result_file << '\n';
#ifdef QUALITY
    if (!settings.log_file.empty()) {
        out << "Log file: " << settings.log_file << '\n';
    }
#endif
#ifdef USE_UINT8
    out << "Data type: uint8_t\n";
#else
    out << "Data type: unsigned long\n";
#endif
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

void generate_keys(Settings const& settings, thread_coordination::Context const& ctx, std::vector<key_type>& keys,
                   std::vector<key_type>& prefill) {
    auto id = ctx.get_id();
    std::seed_seq seed{settings.seed, id};
    std::default_random_engine rng(seed);

    auto first_key = keys.begin() + id * static_cast<std::ptrdiff_t>(settings.elements_per_thread);
    auto last_key = first_key + static_cast<std::ptrdiff_t>(settings.elements_per_thread);
    auto key_dist = std::uniform_int_distribution<key_type>(0, settings.max_key);
    ctx.execute_synchronized(
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

    prefill.resize(settings.prefill_per_thread);
    std::generate(std::begin(prefill), std::end(prefill), [&key_dist, &rng]() { return key_dist(rng); });
}

void prefill(thread_coordination::Context const& ctx, handle_type& handle, std::vector<key_type> const& prefill) {
    ctx.execute_synchronized([&]() {
        for (auto const& key : prefill) {
            handle.push(key);
        }
    });
}

auto execute_mixed(thread_coordination::Context const& ctx, handle_type& handle, std::vector<key_type> const& keys) {
    return ctx.execute_synchronized_blockwise(keys.size(), [&](auto start, auto count) {
        for (auto i = start; i < start + count; ++i) {
            handle.push(keys[i]);
            while (!handle.try_pop()) {
            }
        }
    });
}

auto execute_split_push(thread_coordination::Context const& ctx, handle_type& handle,
                        std::vector<key_type> const& keys) {
    return ctx.execute_synchronized_blockwise(keys.size(), [&](auto start, auto count) {
        for (auto i = start; i < start + count; ++i) {
            handle.push(keys[i]);
        }
    });
}

auto execute_split_pop(thread_coordination::Context const& ctx, handle_type& handle, std::atomic_llong& pop_count) {
    return ctx.execute_synchronized([&]() {
        long long elements_left;
        do {
            long long num_pops = 0;
            while (handle.try_pop()) {
                ++num_pops;
            }
            elements_left = num_pops == 0 ? pop_count.load(std::memory_order_relaxed)
                                          : pop_count.fetch_sub(num_pops, std::memory_order_relaxed) - num_pops;
        } while (elements_left > 0);
    });
}

auto execute_increment(thread_coordination::Context const& ctx, handle_type& handle,
                       std::vector<key_type> const& keys) {
    return ctx.execute_synchronized_blockwise(keys.size(), [&](auto start_index, auto count) {
        for (auto i = start_index; i < start_index + count; ++i) {
            do {
                auto retval = handle.try_pop();
                if (retval) {
                    handle.push(retval->first + keys[i]);
                    break;
                }
            } while (true);
        }
    });
}

struct ThreadData {
    thread_coordination::time_result_type work_time{};
#ifdef WITH_PAPI
    std::vector<long long> papi_counter{};
#endif
#ifdef MQ_COUNT_STATS
    multiqueue::Counters pq_stats{};
#endif
#ifdef QUALITY
    std::vector<OperationLog> op_log{};
#else
    OperationCounters op_counters{};
#endif
};

void write_csv(std::vector<ThreadData> const& stats, std::ostream& out) {
    for (const auto& s : stats) {
        // clang-format off
        out << s.work_time.start.time_since_epoch().count() << ','
            << s.work_time.end.time_since_epoch().count() << ','
            << s.op_counters.push << ','
            << s.op_counters.pop << ','
            << s.op_counters.failed_pop << ',';
        // clang-format on
#ifdef MQ_COUNT_STATS
        // clang-format off
        out << ','
            << s.pq_stats.locked_push_pq << ','
            << s.pq_stats.empty_pop_pqs << ','
            << s.pq_stats.locked_pop_pq << ','
            << s.pq_stats.stale_pop_pq;
        // clang-format on
#else
        out << ",0,0,0,0,0,0,0";
#endif
#ifdef WITH_PAPI
        for (auto c : s.papi_counter) {
            out << ',' << c;
        }
#endif
        out << '\n';
    }
}

ThreadData execute_benchmark(Settings const& settings, thread_coordination::Context const& ctx, handle_type& handle,
                             std::vector<key_type> const& keys) {
    ThreadData thread_data;
#ifdef MQ_COUNT_STATS
    auto pq_stats = handle.get_stats();
#endif
#ifdef WITH_PAPI
    int event_set = PAPI_NULL;
    bool papi_started = false;
    if (!settings.papi_events.empty()) {
        try {
            event_set = start_papi(settings.papi_events);
            papi_started = true;
        } catch (std::exception const& e) {
            ctx.write(std::cerr) << "Error: " << e.what() << '\n';
            std::fill(thread_data.papi_counter.begin(), thread_data.papi_counter.end(), -1);
        }
    }
#endif

    switch (settings.work_mode) {
        case WorkMode::Mixed:
            thread_data.work_time = execute_mixed(ctx, handle, keys);
            break;
        case WorkMode::Split:
            if (ctx.get_id() < settings.num_push_threads) {
                thread_data.work_time = execute_split_push(ctx, handle, keys);
            } else {
                static std::atomic_llong pop_count =
                    static_cast<long long>(settings.prefill_per_thread + settings.elements_per_thread) *
                    settings.num_threads;
                thread_data.work_time = execute_split_pop(ctx, handle, pop_count);
            }
            break;
        case WorkMode::Increment:
            thread_data.work_time = execute_increment(ctx, handle, keys);
            break;
    }

#ifdef WITH_PAPI
    if (papi_started) {
        if (int ret = PAPI_stop(event_set, thread_data.papi_counter.data()); ret != PAPI_OK) {
            ctx.write(std::cerr) << "Error: Failed to stop performance counters\n";
            std::fill(thread_data.papi_counter.begin(), thread_data.papi_counter.end(), -1);
        }
    }
#endif
#ifdef MQ_COUNT_STATS
    thread_data.pq_stats = handle.get_stats() - pq_stats;
#endif
#ifdef QUALITY
    thread_data.op_log = std::move(handle).get_log();
#else
    thread_data.op_counters = handle.get_counters();
#endif
}

std::vector<ThreadData> run_benchmark(Settings const& settings) {
    auto pq =
        pq_type(settings.num_threads, settings.prefill_per_thread * static_cast<std::size_t>(settings.num_threads),
                settings.pq_settings);
    std::cout << "Priority queue: ";
    pq.describe(std::cout) << '\n';
    auto keys = std::vector<key_type>(static_cast<std::size_t>(settings.num_threads) * settings.elements_per_thread);
    std::vector<ThreadData> results(static_cast<std::size_t>(settings.num_threads));
    thread_coordination::TaskHandle task_handle(settings.num_threads, [&](auto ctx) {
        auto handle = Handle(ctx.get_id(), pq);
        if (ctx.get_id() == 0) {
            std::clog << "Generating keys..." << std::flush;
            auto t_start = std::chrono::steady_clock::now();
            generate_keys(settings, ctx, keys);
            auto t_end = std::chrono::steady_clock::now();
            std::clog << "done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count()
                      << "ms)" << std::endl;
            std::clog << "Prefilling... " << std::flush;
            t_start = std::chrono::steady_clock::now();
            prefill(ctx, handle, keys);
            t_end = std::chrono::steady_clock::now();
            std::clog << "done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count()
                      << "ms)" << std::endl;
            std::clog << "Running benchmark..." << std::flush;
            t_start = std::chrono::steady_clock::now();
            results[ctx.get_id()] = execute_benchmark(settings, ctx, handle, keys);
            t_end = std::chrono::steady_clock::now();
            std::clog << "done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count()
                      << "ms)" << std::endl;
        } else {
            generate_keys(settings, ctx, keys);
            prefill(ctx, handle, keys);
            results[ctx.get_id()] = execute_benchmark(settings, ctx, handle, keys);
        }
    });
    task_handle.wait();
    return results;
}

int main(int argc, char* argv[]) {
    std::cout << get_build_info();
    Settings settings;
    cxxopts::Options cmd(argv[0]);
    cmd.add_options()
        // clang-format off
        ("h,help", "Print this help", cxxopts::value<bool>(settings.help))
        ("j,threads", "The number of threads", cxxopts::value<int>(settings.num_threads), "NUMBER")
        ("p,prefill", "The prefill per thread", cxxopts::value<std::size_t>(settings.prefill_per_thread), "NUMBER")
        ("n,keys", "The number of keys per thread", cxxopts::value<std::size_t>(settings.elements_per_thread), "NUMBER")
        ("w,work-mode", "Specify the work mode ([m]ixed, [s]plit, [i]ncrement)", cxxopts::value<char>(settings.work_mode), "STRING")
        ("i,push-threads", "The number of pushing threads in split mode", cxxopts::value<int>(settings.num_push_threads), "NUMBER")
        ("e,element-distribution", "Specify the element distribution ([u]niform, [a]scending, [d]escending)", cxxopts::value<char>(settings.key_distribution), "STRING")
        ("m,max", "Specify the max key", cxxopts::value<key_type>(settings.max_key), "NUMBER")
        ("s,seed", "Specify the initial seed", cxxopts::value<int>(settings.seed), "NUMBER")
        ("x,result-file", "Path to write per-thread results to in CSV format", cxxopts::value<std::filesystem::path>(settings.result_file), "PATH")
#ifdef QUALITY
        ("histogram-file", "Path to write the histogram to", cxxopts::value<std::filesystem::path>(settings.histogram_file), "PATH")
        ("l,log-file", "Path to write the operation log to", cxxopts::value<std::filesystem::path>(settings.log_file), "PATH")
#endif
#ifdef WITH_PAPI
        ("r,pc", "Enable performance counter", cxxopts::value<>(settings.enable_performance_counter))
#endif
        // clang-format on
        ;
    pq_type::add_options(cmd, settings.pq_settings);
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
    if (!settings->papi_events.empty()) {
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

    auto result_file = std::ofstream(settings->result_file);
    if (!result_file) {
        std::cerr << "Error: Could not open result file " << settings->result_file << " for writing" << std::endl;
        return EXIT_FAILURE;
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
#endif

    auto results = run_benchmark(*settings);

#ifdef QUALITY
    if (!settings->result_file.empty()) {
        write_result_csv(results, result_file);
    }

#ifdef QUALITY
    std::vector<operation_log::OperationLog> logs;
    for (auto& result : results) {
        logs.push_back(std::move(log));
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
