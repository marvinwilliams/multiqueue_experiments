#include "build_info.hpp"
#include "priority_queue_selection.hpp"
#include "task.hpp"

#include "cxxopts.hpp"

#ifdef QUALITY
#undef USE_UINT8
#include "operation_log.hpp"
#else
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
static constexpr auto max_key = std::numeric_limits<key_type>::max();
#else
using key_type = unsigned long;
using value_type = unsigned long;
static constexpr auto max_key = static_cast<key_type>(1) << 30;
#endif

using pq_type = PriorityQueue<key_type, value_type, true>;
using handle_type = pq_type::handle_type;

template <typename Integer = long long>
class ChunkExecutor {
    std::atomic<Integer> counter_;
    Integer chunk_size_;
    Integer max_;

   public:
    explicit ChunkExecutor(Integer chunk_size, Integer max) : chunk_size_{chunk_size}, max_{max} {
    }

    template <typename Work, typename... Args>
    std::chrono::nanoseconds operator()(task::Control const& tc, Work work, Args&&... args) {
        std::chrono::nanoseconds total_time{};
        if (tc.id() == 0) {
            counter_.store(0, std::memory_order_relaxed);
        }
        tc.synchronize();
        auto from = counter_.fetch_add(chunk_size_, std::memory_order_relaxed);
        while (from < max_) {
            auto to = std::min(from + chunk_size_, max_);
            auto t_start = std::chrono::high_resolution_clock::now();
            for (; from < to; ++from) {
                work(from, args...);  // no perfect forwarding here
            }
            auto t_end = std::chrono::high_resolution_clock::now();
            total_time += t_end - t_start;
            from = counter_.fetch_add(chunk_size_, std::memory_order_relaxed);
        }
        tc.synchronize();
        return total_time;
    }
};

template <typename Integer = long long>
class SingleExecutor {
    std::atomic<Integer> counter_;
    Integer max_;

   public:
    explicit SingleExecutor(Integer max) : max_{max} {
    }

    template <typename Work, typename... Args>
    std::chrono::nanoseconds operator()(task::Control const& tc, Work work, Args&&... args) {
        std::chrono::nanoseconds total_time{};
        if (tc.id() == 0) {
            counter_.store(0, std::memory_order_relaxed);
        }
        tc.synchronize();
        for (Integer i = counter_.fetch_add(1, std::memory_order_relaxed); i < max_;
             i = counter_.fetch_add(1, std::memory_order_relaxed)) {
            auto t_start = std::chrono::high_resolution_clock::now();
            work(i, args...);  // no perfect forwarding here
            auto t_end = std::chrono::high_resolution_clock::now();
            total_time += t_end - t_start;
        }
        tc.synchronize();
        return total_time;
    }
};

template <typename Integer = long long>
class SplitExecutor {
    Integer max_;

   public:
    explicit SplitExecutor(Integer max) : max_{max} {
    }

    template <typename Work, typename... Args>
    std::chrono::nanoseconds operator()(task::Control const& tc, Work work, Args&&... args) {
        auto from = (static_cast<Integer>(tc.id()) * max_) / static_cast<Integer>(tc.num_threads());
        auto to = tc.id() == tc.num_threads() - 1 ? max_ : from + (max_ / static_cast<Integer>(tc.num_threads()));
        tc.synchronize();
        auto t_start = std::chrono::high_resolution_clock::now();
        for (; from != to; ++from) {
            work(from, args...);  // no perfect forwarding here
        }
        auto t_end = std::chrono::high_resolution_clock::now();
        tc.synchronize();
        return t_end - t_start;
    }
};

class Benchmark {
    enum class Mode { Increment, Alternate, Push, Pop };
    enum class KeyDistribution { Random, Ascending, Descending };

    struct Settings {
        int num_threads = 4;
        long long prefill_per_thread = 1 << 20;
        long long elements_per_thread = 1 << 24;
        Mode mode = Mode::Increment;
        KeyDistribution key_distribution = KeyDistribution::Random;
        pq_type::config_type pq_settings{};
        key_type min_prefill = 1;
        key_type max_prefill = max_key;
        key_type min_increment = 1;
        key_type max_increment = 100;
        key_type min_random = 1;
        key_type max_random = max_key;
        long long block_size = 1 << 12;
        int seed = 1;
#ifdef QUALITY
        std::filesystem::path log_file;
#endif
#ifdef USE_PAPI
        std::vector<std::string> papi_events;
#endif
    };

    struct BenchmarkStats {
        std::chrono::nanoseconds time{};
#ifdef USE_PAPI
        std::vector<long long> papi_counter{};
#endif
    };

    struct OperationCount {
        long long pop = 0;
        long long push = 0;
        long long failed_pop = 0;
    };

    struct Stats {
        BenchmarkStats benchmark_stats{};
        OperationCount op_count{};
#ifdef MQ_COUNT_STATS
        multiqueue::Counters mq_stats{};
#endif
    };

   private:
    Settings settings_;
    std::vector<key_type> keys_;

#ifdef USE_PAPI
    int start_papi() const {
        int event_set = PAPI_NULL;
        if (int ret = PAPI_register_thread(); ret != PAPI_OK) {
            throw std::runtime_error("Failed to register thread for PAPI");
        }
        if (int ret = PAPI_create_eventset(&event_set); ret != PAPI_OK) {
            throw std::runtime_error("Failed to create PAPI event set");
        }
        for (auto const& name : settings_.papi_events) {
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

    void generate_keys(int id) {
        std::seed_seq seed{settings_.seed, id};
        std::default_random_engine rng(seed);

        auto prefill_begin = keys_.begin() + id * settings_.prefill_per_thread;
        auto prefill_end = prefill_begin + settings_.prefill_per_thread;
        auto prefill_dist = std::uniform_int_distribution<key_type>(settings_.min_prefill, settings_.max_prefill);
        std::generate(prefill_begin, prefill_end, [&prefill_dist, &rng]() { return prefill_dist(rng); });
        auto workload_begin =
            keys_.begin() + settings_.num_threads * settings_.prefill_per_thread + id * settings_.elements_per_thread;
        auto workload_end = workload_begin + settings_.elements_per_thread;
        switch (settings_.mode) {
            case Mode::Increment: {
                auto dist = std::uniform_int_distribution<key_type>(settings_.min_increment, settings_.max_increment);
                std::generate(workload_begin, workload_end, [&dist, &rng]() { return dist(rng); });
                break;
            }
            case Mode::Push:
            case Mode::Alternate: {
                switch (settings_.key_distribution) {
                    case KeyDistribution::Random: {
                        auto dist = std::uniform_int_distribution<key_type>(settings_.min_random, settings_.max_random);
                        std::generate(workload_begin, workload_end, [&dist, &rng]() { return dist(rng); });
                        break;
                    }
                    case KeyDistribution::Ascending: {
                        std::iota(workload_begin, workload_end, std::distance(keys_.begin(), workload_begin));
                        break;
                    }
                    case KeyDistribution::Descending: {
                        std::iota(std::reverse_iterator(workload_end), std::reverse_iterator(workload_begin),
                                  std::distance(workload_end, keys_.end()));
                        break;
                    }
                }
            }
            case Mode::Pop:
                break;
        }
    }

    void prefill(task::Control const& tc, handle_type& h) {
        auto offset = static_cast<std::size_t>(tc.id() * settings_.prefill_per_thread);
        for (std::size_t i = 0; i < static_cast<std::size_t>(settings_.prefill_per_thread); ++i) {
#ifdef QUALITY
            h.push({keys_[offset + i], offset + i});
#else
            h.push({keys_[offset + i], keys_[offset + i]});
#endif
        }
    }

    void write_stats(Stats const& stats, std::ostream& out) {
        out << stats.benchmark_stats.time.count() << ',';
#ifdef USE_PAPI
        for (auto c : stats.benchmark_stats.papi_counter) {
            out << ',' << c;
        }
#endif
        // clang-format off
        out << stats.op_count.push << ','
            << stats.op_count.pop << ','
            << stats.op_count.failed_pop;
        // clang-format on
#ifdef MQ_COUNT_STATS
        // clang-format off
        out << ','
            << stats.mq_stats.locked_push_pq << ','
            << stats.mq_stats.empty_pop_pqs << ','
            << stats.mq_stats.locked_pop_pq << ','
            << stats.mq_stats.stale_pop_pq;
        // clang-format on
#endif
        out << '\n';
    }

    template <typename Executor>
    Stats work(task::Control const& tc, handle_type& h, Executor& executor) {
        OperationCount op_count{};
#ifdef QUALITY
        auto push = [&h, &op_count, this](key_type key, std::size_t index) {
            h.push(key, index);
            auto tick = get_tick();
            push_log[index] = {tick, key};
            ++op_count.push;
        };

        auto pop = [&h, &op_count, this](std::size_t index) {
            auto tick = get_tick();
            auto retval = h.try_pop();
            while (!retval) {
                ++op_count.failed_pop;
                tick = get_tick();
                retval = h.try_pop();
            }
            pop_log[index] = {tick, retval->second};
            ++op_count.pop;
            return retval->first;
        };
#else
        auto push = [&h, &op_count](auto key, std::size_t) {
            h.push((key), (key));
            ++op_count.push;
        };

        auto pop = [&h, &op_count](std::size_t) {
            auto retval = h.try_pop();
            while (!retval) {
                ++op_count.failed_pop;
                retval = h.try_pop();
            }
            ++op_count.pop;
            return retval->first;
        };
#endif
        auto offset = static_cast<std::size_t>(settings_.num_threads * settings_.prefill_per_thread);
        Stats stats{};
#ifdef MQ_COUNT_STATS
        h.reset_counters();
#endif
        switch (settings_.mode) {
            case Mode::Alternate: {
                auto work = [offset, this](std::size_t index) {
                    push(keys_[offset + index], offset + index);
                    pop(index);
                };
                stats.benchmark_stats = execute(tc, executor, work);
                break;
            }
            case Mode::Increment: {
                auto work = [offset, this](std::size_t index) {
                    auto key = pop(index);
                    push(key + keys_[offset + index], offset + index);
                };
                stats.benchmark_stats = execute(tc, executor, work);
                break;
            }
            case Mode::Push: {
                auto work = [offset, this](std::size_t index) { push(keys_[offset + index], offset + index); };
                stats.benchmark_stats = execute(tc, executor, work);
                break;
            }
            case Mode::Pop: {
                auto work = [this](std::size_t index) { pop(index); };
                stats.benchmark_stats = execute(tc, executor, work);
                break;
            }
        }
        stats.op_count = op_count;
#ifdef MQ_COUNT_STATS
        stats.mq_stats = h.get_counters();
#endif
        return stats;
    }

    template <typename Executor, typename Work>
    BenchmarkStats execute(task::Control const& tc, Executor& executor, Work work) {
        BenchmarkStats stats{};
#ifdef USE_PAPI
        int event_set = PAPI_NULL;
        bool papi_started = false;
        if (!settings_.papi_events.empty()) {
            stats.papi_counter.resize(settings_.papi_events.size());
            try {
                event_set = start_papi();
                papi_started = true;
            } catch (std::exception const& e) {
                tc.write(std::cerr) << "Error: " << e.what() << '\n';
                std::fill(stats.papi_counter.begin(), stats.papi_counter.end(), -1);
            }
        }
#endif
        stats.time = executor(tc, work);

#ifdef USE_PAPI
        if (papi_started) {
            if (int ret = PAPI_stop(event_set, stats.papi_counter.data()); ret != PAPI_OK) {
                tc.write(std::cerr) << "Error: Failed to stop performance counters\n";
                std::fill(stats.papi_counter.begin(), stats.papi_counter.end(), -1);
            }
        }
#endif
        return stats;
    }

   public:
    Benchmark(Settings settings)
        : settings_(settings),
          keys_(static_cast<std::size_t>(settings_.num_threads *
                                         (settings_.prefill_per_thread + settings_.elements_per_thread))) {
    }

    template <typename Executor>
    Stats benchmark_thread(task::Control const& tc, handle_type h, Executor& executor) {
        Stats stats{};
        std::chrono::steady_clock::time_point t_start;
        std::chrono::steady_clock::time_point t_end;
        tc.once([]() { std::clog << "Generating keys..." << std::flush; });
        tc.synchronize();
        tc.once([&t_start]() { t_start = std::chrono::steady_clock::now(); });
        tc.synchronize();
        generate_keys(tc.id());
        tc.synchronize();
        tc.once([&t_end, &t_start]() {
            t_end = std::chrono::steady_clock::now();
            std::clog << "done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count()
                      << "ms)" << std::endl;
            std::clog << "Prefilling..." << std::flush;
            t_start = std::chrono::steady_clock::now();
        });
        tc.synchronize();
        prefill(tc, h);
        tc.synchronize();
        tc.once([&t_end, &t_start]() {
            t_end = std::chrono::steady_clock::now();
            std::clog << "done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count()
                      << "ms)" << std::endl;
            std::clog << "Running benchmark..." << std::flush;
            t_start = std::chrono::steady_clock::now();
        });
        stats = work(tc, h, executor);
        tc.once([&t_end, &t_start]() {
            t_end = std::chrono::steady_clock::now();
            std::clog << "done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count()
                      << "ms)" << std::endl;
        });
        return stats;
    }

    void run() {
        auto pq = pq_type(settings_.num_threads,
                          static_cast<std::size_t>(settings_.prefill_per_thread * settings_.num_threads),
                          settings_.pq_settings);
        std::clog << "Priority queue: ";
        pq.describe(std::clog) << '\n';

        std::vector<Stats> stats(static_cast<std::size_t>(settings_.num_threads));
        affinity::NUMA numa_affinity{CORES_PER_NUMA_NODE, NUM_NUMA_NODES};
    // TODO: Executor
        task::Runner runner(numa_affinity, settings_.num_threads, [&](task::Control tc) {
            auto h = pq.get_handle();
            benchmark_thread();
        });
        task_handle.wait();
        return results;
    }

    BenchmarkMode parse_work_mode(char c) {
        switch (c) {
            case 'm':
                return BenchmarkMode::Mixed;
            case 'u':
                return BenchmarkMode::Push;
            case 'o':
                return BenchmarkMode::Pop;
            case 'i':
                return BenchmarkMode::Increment;
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
};

bool validate_settings(Benchmark::Settings const& settings) if (settings.num_threads <= 0) {
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

[[nodiscard]] auto mode_str(Benchmark::Mode mode) {
    switch (mode) {
        case Mode::Alternate:
            return "alternate";
        case Mode::Increment:
            return "increment";
        case Mode::Push:
            return "push";
        case Mode::Pop:
            return "pop";
    }
    return "unknown";
}
[[nodiscard]] auto key_dist_str(Benchmark::KeyDistribution dist) {
    switch (dist) {
        case Benchmark::KeyDistribution::Uniform:
            return "uniform";
        case Benchmark::KeyDistribution::Ascending:
            return "ascending";
        case Benchmark::KeyDistribution::Descending:
            return "descending";
    }
    return "unknown";
}

void write_settings(Benchmark::Settings const& settings, std::ostream& out) {
    out << "Threads: " << settings_.num_threads << '\n'
        << "Prefill per thread: " << settings_.prefill_per_thread << '\n'
        << "Elements per thread: " << settings_.elements_per_thread << '\n'
        << "Mode: " << work_mode_name(settings_.mode) << '\n'
        << "Max key: " << static_cast<unsigned long>(settings_.max_key) << '\n'
        << "Seed: " << settings_.seed << '\n';
#ifdef USE_PAPI
    if (!settings_.papi_events.empty()) {
        out << "PAPI events: ";
        std::copy(settings_.papi_events.begin(), settings_.papi_events.end(),
                  std::ostream_iterator<std::string>(out, " "));
        out << '\n';
    }
#endif
#ifdef QUALITY
    if (!settings_.log_file.empty()) {
        out << "Log operations to: " << settings_.log_file << '\n';
    }
#else
#ifdef USE_UINT8
    out << "Data type: uint8_t\n";
#else
    out << "Data type: unsigned long\n";
#endif
#endif
}
int main(int argc, char* argv[]) {
    write_build_info(std::clog);
    Settings settings;
    cxxopts::Options cmd(argv[0]);
    cmd.add_options()
        // clang-format off
        ("h,help", "Print this help", cxxopts::value<bool>())
        ("j,threads", "The number of threads", cxxopts::value<int>(settings_.num_threads), "NUMBER")
        ("p,prefill", "The prefill per thread", cxxopts::value<std::size_t>(settings_.prefill_per_thread), "NUMBER")
        ("n,keys", "The number of keys per thread", cxxopts::value<std::size_t>(settings_.elements_per_thread), "NUMBER")
        ("w,work-mode", "Specify the work mode ([m]ixed, p[u]sh, p[o]p, [i]ncrement)", cxxopts::value<char>(), "STRING")
        ("d,key-dist", "Specify the key distribution ([u]niform, [a]scending, [d]escending)", cxxopts::value<char>(), "STRING")
        ("m,max", "Specify the max key", cxxopts::value<key_type>(settings_.max_key), "NUMBER")
        ("s,seed", "Specify the initial seed", cxxopts::value<int>(settings_.seed), "NUMBER")
#ifdef QUALITY
        ("l,log-file", "Path to write the operation log to", cxxopts::value<std::filesystem::path>(settings_.log_file), "PATH")
#endif
#ifdef USE_PAPI
        ("r,pc", "Performance counters", cxxopts::value<std::vector<std::string>>(settings_.papi_events))
#endif
        // clang-format on
        ;
    pq_type::add_options(cmd, settings_.pq_settings);
    try {
        auto args = cmd.parse(argc, argv);
        if (args.count("help") > 0) {
            std::clog << cmd.help() << std::endl;
            return EXIT_SUCCESS;
        }
        if (args.count("work-mode") > 0) {
            settings_.work_mode = parse_work_mode(args["work-mode"].as<char>());
        }
        if (args.count("key-dist") > 0) {
            settings_.key_distribution = parse_key_distribution(args["key-dist"].as<char>());
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
    if (!settings_.papi_events.empty()) {
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
    if (!settings_.log_file.empty()) {
        log_file = std::ofstream(settings_.log_file);
        if (!log_file) {
            std::cerr << "Error: Could not open log file " << settings_.log_file << " for writing" << std::endl;
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
    if (!settings_.log_file.empty()) {
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
    out << "thread,start_time,end_time,pushes,pops,failed_pops,locked_push_pq,empty_pop_pqs,"
           "locked_pop_pq,stale_pop_pq";
#ifdef USE_PAPI
    for (auto const& e : settings_.papi_events) {
        out << ',' << e;
    }
#endif
    out << '\n';
    write_thread_data(results, std::cout);
#endif
}
