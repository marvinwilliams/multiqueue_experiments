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

#ifdef CORES_PER_NUMA_NODE
static constexpr auto cores_per_numa_node = CORES_PER_NUMA_NODE;
#else
static constexpr auto cores_per_numa_node = 4;
#endif
#ifdef NUM_NUMA_NODES
static constexpr auto num_numa_nodes = NUM_NUMA_NODES;
#else
static constexpr auto num_numa_nodes = 16;
#endif

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

template <typename Integer = std::size_t>
class ChunkExecutor {
    std::atomic<Integer> counter_;
    Integer chunk_size_;

   public:
    explicit ChunkExecutor(Integer chunk_size) : chunk_size_{chunk_size} {
    }

    template <typename Work, typename... Args>
    auto operator()(task::Control const& tc, Integer max, Work work, Args&&... args) {
        if (tc.id() == 0) {
            counter_.store(0, std::memory_order_relaxed);
        }
        tc.synchronize();
        auto t_start = std::chrono::high_resolution_clock::now();
        auto from = counter_.fetch_add(chunk_size_, std::memory_order_relaxed);
        while (from < max) {
            auto to = std::min(from + chunk_size_, max);
            for (; from < to; ++from) {
                work(from, args...);  // no perfect forwarding here
            }
            from = counter_.fetch_add(chunk_size_, std::memory_order_relaxed);
        }
        auto t_end = std::chrono::high_resolution_clock::now();
        tc.synchronize();
        return std::pair{t_start, t_end};
    }
};

template <typename Integer = std::size_t>
class SingleExecutor {
    std::atomic<Integer> counter_;

   public:
    template <typename Work, typename... Args>
    auto operator()(task::Control const& tc, Integer max, Work work, Args&&... args) {
        if (tc.id() == 0) {
            counter_.store(0, std::memory_order_relaxed);
        }
        tc.synchronize();
        auto t_start = std::chrono::high_resolution_clock::now();
        for (Integer i = counter_.fetch_add(1, std::memory_order_relaxed); i < max;
             i = counter_.fetch_add(1, std::memory_order_relaxed)) {
            work(i, args...);  // no perfect forwarding here
        }
        auto t_end = std::chrono::high_resolution_clock::now();
        tc.synchronize();
        return std::pair{t_start, t_end};
    }
};

template <typename Integer = std::size_t>
class SplitExecutor {
   public:
    template <typename Work, typename... Args>
    auto operator()(task::Control const& tc, Integer max, Work work, Args&&... args) {
        auto from = (static_cast<Integer>(tc.id()) * max) / static_cast<Integer>(tc.num_threads());
        auto to = tc.id() == tc.num_threads() - 1 ? max : from + (max / static_cast<Integer>(tc.num_threads()));
        tc.synchronize();
        auto t_start = std::chrono::high_resolution_clock::now();
        for (; from != to; ++from) {
            work(from, args...);  // no perfect forwarding here
        }
        auto t_end = std::chrono::high_resolution_clock::now();
        tc.synchronize();
        return std::pair{t_start, t_end};
    }
};

class Benchmark {
   public:
    enum class Mode { Increment, Alternate, Push, Pop };
    enum class KeyDistribution { Random, Ascending, Descending };
    enum class ExecutionMode { Chunk, Single, Split };

    struct Settings {
        int num_threads = 4;
        long long prefill_per_thread = 1 << 20;
        long long operations_per_thread = 1 << 24;
        Mode mode = Mode::Increment;
        KeyDistribution key_distribution = KeyDistribution::Random;
        ExecutionMode execution_mode = ExecutionMode::Chunk;
        pq_type::config_type pq_settings{};
        key_type min_prefill = 1;
        key_type max_prefill = max_key;
        key_type min_increment = 1;
        key_type max_increment = 100;
        key_type min_random = 1;
        key_type max_random = max_key;
        long long chunk_size = 1 << 12;
        int seed = 1;
#ifdef QUALITY
        std::filesystem::path log_file;
#endif
#ifdef USE_PAPI
        std::vector<std::string> papi_events;
#endif
    };

    class ValidSettings {
        Settings settings_{};

       public:
        explicit ValidSettings(Settings settings) : settings_{std::move(settings)} {
            if (settings.num_threads <= 0) {
                throw std::runtime_error{"Error: Number of threads must be greater than 0"};
            }
            if (settings.min_prefill <= 0 || settings.min_random <= 0) {
                throw std::runtime_error{"Error: Keys must be greater than 0"};
            }
            if (settings.max_prefill < settings.min_prefill || settings.max_random < settings.min_random ||
                settings.max_increment < settings.min_increment) {
                throw std::runtime_error{"Error: Max must be greater than min"};
            }
            if (settings.chunk_size <= 0) {
                throw std::runtime_error{"Error: Chunk size must be greater than 0"};
            }
            if (settings.prefill_per_thread < 0 || settings.operations_per_thread < 0) {
                throw std::runtime_error{"Error: Number of operations must be greater than 0"};
            }
#ifdef USE_PAPI
            for (auto const& name : settings.papi_events) {
                if (PAPI_query_named_event(name.c_str()) != PAPI_OK) {
                    throw std::runtime_error{"Error: Invalid PAPI event name: " + name};
                }
            }
#endif
        }

        explicit ValidSettings() : ValidSettings{Settings{}} {
        }
        [[nodiscard]] Settings const& get() const& {
            return settings_;
        }
        Settings&& get() && {
            return std::move(settings_);
        }
    };

    struct BenchmarkStats {
        std::chrono::high_resolution_clock::time_point start_time{};
        std::chrono::high_resolution_clock::time_point end_time{};
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
    std::vector<Stats> stats_;
#ifdef QUALITY
    operation_log::OperationLog<key_type> op_log_;
#endif

    static int pushes_per_iteration(Mode mode) {
        switch (mode) {
            case Mode::Alternate:
            case Mode::Increment:
            case Mode::Push:
                return 1;
            case Mode::Pop:
                return 0;
        }
    }

    static int pops_per_iteration(Mode mode) {
        switch (mode) {
            case Mode::Alternate:
            case Mode::Increment:
            case Mode::Pop:
                return 1;
            case Mode::Push:
                return 0;
        }
    }

    static long long num_iterations(Settings const& settings) {
        return settings.num_threads * settings.operations_per_thread /
            (pushes_per_iteration(settings.mode) + pops_per_iteration(settings.mode));
    }

    static long long num_pushes(Settings const& settings) {
        return settings.num_threads * settings.prefill_per_thread +
            pushes_per_iteration(settings.mode) * num_iterations(settings);
    }

    static long long num_pops(Settings const& settings) {
        return pops_per_iteration(settings.mode) * num_iterations(settings);
    }

    static long long max_num_elements(Settings const& settings) {
        return settings.num_threads *
            (settings.prefill_per_thread + (settings.mode == Mode::Push ? settings.operations_per_thread : 0));
    }

#ifdef USE_PAPI
    [[nodiscard]] int start_papi() const {
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
        auto first_workload = keys_.begin() + settings_.num_threads * settings_.prefill_per_thread;
        auto workload_begin =
            first_workload + id * (std::distance(first_workload, keys_.end()) / settings_.num_threads);
        auto workload_end =
            first_workload + (id + 1) * (std::distance(first_workload, keys_.end()) / settings_.num_threads);
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
            op_log_.pushes[offset + i] = {0, keys_[offset + i]};
#else
            h.push({keys_[offset + i], keys_[offset + i]});
#endif
        }
    }

    template <typename Executor>
    void work(task::Control const& tc, handle_type& h, Executor& executor) {
        OperationCount op_count{};
#ifdef QUALITY
        auto push = [&h, &op_count, this](key_type key, std::size_t index) {
            h.push(key, static_cast<value_type>(index));
            auto tick = operation_log::get_tick();
            op_log_.pushes[index] = {tick, key};
            ++op_count.push;
        };

        auto pop = [&h, &op_count, this](std::size_t index) {
            auto tick = operation_log::get_tick();
            auto retval = h.try_pop();
            while (!retval) {
                ++op_count.failed_pop;
                tick = operation_log::get_tick();
                retval = h.try_pop();
            }
            op_log_.pop_log[index] = {tick, static_cast<std::size_t>(retval->second)};
            ++op_count.pop;
            return retval->first;
        };
#else
        auto push = [&h, &op_count](auto key, std::size_t) {
            h.push({key, key});
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
#ifdef MQ_COUNT_STATS
        h.reset_counters();
#endif
        switch (settings_.mode) {
            case Mode::Alternate: {
                auto work = [&, offset](std::size_t index) {
                    push(keys_[offset + index], offset + index);
                    pop(index);
                };
                execute(tc, executor, work);
                break;
            }
            case Mode::Increment: {
                auto work = [&, offset](std::size_t index) {
                    auto key = pop(index);
                    push(key + keys_[offset + index], offset + index);
                };
                execute(tc, executor, work);
                break;
            }
            case Mode::Push: {
                auto work = [&, offset](std::size_t index) { push(keys_[offset + index], offset + index); };
                execute(tc, executor, work);
                break;
            }
            case Mode::Pop: {
                auto work = [&](std::size_t index) { pop(index); };
                execute(tc, executor, work);
                break;
            }
        }
        stats_[static_cast<std::size_t>(tc.id())].op_count = op_count;
#ifdef MQ_COUNT_STATS
        stats_[static_cast<std::size_t>(tc.id())].mq_stats = h.get_counters();
#endif
    }

    template <typename Executor, typename Work>
    void execute(task::Control const& tc, Executor& executor, Work work) {
#ifdef USE_PAPI
        int event_set = PAPI_NULL;
        bool papi_started = false;
        if (!settings_.papi_events.empty()) {
            stats_[tc.id()].papi_counter.resize(settings_.papi_events.size());
            try {
                event_set = start_papi();
                papi_started = true;
            } catch (std::exception const& e) {
                tc.write(std::cerr) << "Error: " << e.what() << '\n';
            }
        }
#endif
        auto iterations =
            static_cast<std::size_t>(settings_.operations_per_thread /
                                     (pushes_per_iteration(settings_.mode) + pops_per_iteration(settings_.mode)));
        auto [t_start, t_end] = executor(tc, iterations, work);

#ifdef USE_PAPI
        if (papi_started) {
            if (int ret = PAPI_stop(event_set, stats_[tc.id()].papi_counter.data()); ret != PAPI_OK) {
                tc.write(std::cerr) << "Error: Failed to stop performance counters\n";
                std::fill(stats_[tc.id()].papi_counter.begin(), stats_[tc.id()].papi_counter.end(), 0);
            }
        }
#endif
        stats_[tc.id()].start_time = t_start;
        stats_[tc.id()].end_time = t_end;
    }

    template <typename Executor>
    void run_thread(task::Control const& tc, pq_type& pq, Executor& executor) {
        auto h = pq.get_handle();
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
        work(tc, h, executor);
        tc.once([&t_end, &t_start]() {
            t_end = std::chrono::steady_clock::now();
            std::clog << "done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count()
                      << "ms)" << std::endl;
        });
    }

    static void write_stats_single(Benchmark::Stats const& stats, std::ostream& out) {
        out << (stats.benchmark_stats.end_time - stats.benchmark_stats.start_time).count() << ',';
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

    void write_stats(std::vector<Benchmark::Stats> const& stats, std::ostream& out) {
        out << "stat,time";
#ifdef USE_PAPI
        for (auto const& name : settings_.papi_events) {
            out << ',' << name;
        }
#endif
        out << ",pushes,pops,failed-pops";
#ifdef MQ_COUNT_STATS
        out << ",locked-push-pq,empty-pop-pqs,locked-pop-pq,stale-pop-pq";
#endif
        out << '\n';
        for (auto const& s : stats) {
            std::cout << "thread";
            write_stats_single(s, out);
        }
        Benchmark::Stats total_stats;
        total_stats.benchmark_stats.start_time = std::chrono::high_resolution_clock::time_point::max();
        total_stats.benchmark_stats.end_time = std::chrono::high_resolution_clock::time_point::min();
        for (auto const& s : stats) {
            total_stats.benchmark_stats.start_time =
                std::min(total_stats.benchmark_stats.start_time, s.benchmark_stats.start_time);
            total_stats.benchmark_stats.end_time =
                std::max(total_stats.benchmark_stats.end_time, s.benchmark_stats.end_time);
#ifdef USE_PAPI
            if (total_stats.benchmark_stats.papi_counter.empty()) {
                total_stats.benchmark_stats.papi_counter = s.benchmark_stats.papi_counter;
            } else {
                std::transform(total_stats.benchmark_stats.papi_counter.begin(),
                               total_stats.benchmark_stats.papi_counter.end(), s.benchmark_stats.papi_counter.begin(),
                               total_stats.benchmark_stats.papi_counter.begin(), std::plus<>());
            }
#endif
            total_stats.op_count.push += s.op_count.push;
            total_stats.op_count.pop += s.op_count.pop;
            total_stats.op_count.failed_pop += s.op_count.failed_pop;
#ifdef MQ_COUNT_STATS
            total_stats.mq_stats.locked_push_pq += s.mq_stats.locked_push_pq;
            total_stats.mq_stats.empty_pop_pqs += s.mq_stats.empty_pop_pqs;
            total_stats.mq_stats.locked_pop_pq += s.mq_stats.locked_pop_pq;
            total_stats.mq_stats.stale_pop_pq += s.mq_stats.stale_pop_pq;
#endif
        }
        out << "total,";
        write_stats_single(total_stats, out);
    }

   public:
    Benchmark(ValidSettings settings)
        : settings_(std::move(settings).get()),
          keys_(static_cast<std::size_t>(num_pushes(settings_))),
          stats_(static_cast<std::size_t>(settings_.num_threads)) {
#ifdef QUALITY
        op_log_.pushes.resize(num_pushes(settings_));
        op_log_.pops.resize(num_pops(settings_));
#endif
    }

    bool run() {
        auto pq =
            pq_type(settings_.num_threads, static_cast<std::size_t>(num_pushes(settings_)), settings_.pq_settings);
        std::clog << "Priority queue: ";
        pq.describe(std::clog) << '\n';
        affinity::NUMA numa_affinity{cores_per_numa_node, num_numa_nodes};
        switch (settings_.execution_mode) {
            case ExecutionMode::Single: {
                SingleExecutor executor{};
                task::Runner runner(numa_affinity, settings_.num_threads, run_thread, std::ref(pq), std::ref(executor));
                break;
            }
            case ExecutionMode::Chunk: {
                ChunkExecutor executor{static_cast<std::size_t>(settings_.chunk_size)};
                task::Runner runner(numa_affinity, settings_.num_threads, run_thread, std::ref(pq), std::ref(executor));
                break;
            }
            case ExecutionMode::Split: {
                SplitExecutor executor;
                task::Runner runner(numa_affinity, settings_.num_threads, run_thread, std::ref(pq), std::ref(executor));
                break;
            }
        }
#ifdef QUALITY
        write_stats(settings_, result.stats, stats_out);
        stats_out.close();

        if (!settings.log_file.empty()) {
            std::clog << "Writing operation logs..." << std::flush;
            auto t_start = std::chrono::steady_clock::now();
            operation_log::write_logs(benchmark.get_operation_logs(), log_out);
            log_out.close();
            auto t_end = std::chrono::steady_clock::now();
            std::clog << "done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count()
                      << "ms)\n";
        }
        if (!operation_log::verify_logs(result.op_log)) {
            std::cerr << "Error: Operation logs are invalid\n";
            return false;
        }
        std::clog << "Replaying logs..." << std::flush;
        auto t_start = std::chrono::steady_clock::now();
        auto metrics = operation_log::replay_logs(result.op_log);
        auto t_end = std::chrono::steady_clock::now();
        std::clog << "done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count()
                  << "ms)\n";
        operation_log::write_metrics(metrics, std::cout);
#else
        write_stats(std::cout);
#endif
    }
};

Benchmark::Mode parse_mode(char c) {
    switch (c) {
        case 'a':
            return Benchmark::Mode::Alternate;
        case 'i':
            return Benchmark::Mode::Increment;
        case 'o':
            return Benchmark::Mode::Pop;
        case 'u':
            return Benchmark::Mode::Push;
        default:
            throw std::invalid_argument("Invalid work mode");
    }
}

Benchmark::KeyDistribution parse_key_distribution(char c) {
    switch (c) {
        case 'r':
            return Benchmark::KeyDistribution::Random;
        case 'a':
            return Benchmark::KeyDistribution::Ascending;
        case 'd':
            return Benchmark::KeyDistribution::Descending;
        default:
            throw std::invalid_argument("Invalid key distribution");
    }
}

Benchmark::ExecutionMode parse_execution_mode(char c) {
    switch (c) {
        case 'c':
            return Benchmark::ExecutionMode::Chunk;
        case 's':
            return Benchmark::ExecutionMode::Single;
        case 'p':
            return Benchmark::ExecutionMode::Split;
        default:
            throw std::invalid_argument("Invalid execution mode");
    }
}

bool validate_settings(Benchmark::Settings const& settings) {
    if (settings.num_threads <= 0) {
        std::cerr << "Error: Number of threads must be greater than 0\n";
        return false;
    }
    if (settings.min_prefill <= 0 || settings.min_random <= 0) {
        std::cerr << "Error: Keys must be greater than 0\n";
        return false;
    }
    if (settings.max_prefill < settings.min_prefill || settings.max_random < settings.min_random ||
        settings.max_increment < settings.min_increment) {
        std::cerr << "Error: Key ranges must be non-empty\n";
        return false;
    }
    if (settings.chunk_size <= 0) {
        std::cerr << "Error: Chunk size must be greater than 0\n";
        return false;
    }
    if (settings.prefill_per_thread < 0 || settings.operations_per_thread < 0) {
        std::cerr << "Error: Prefill and operations must be non-negative\n";
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

[[nodiscard]] auto mode_name(Benchmark::Mode mode) {
    switch (mode) {
        case Benchmark::Mode::Alternate:
            return "alternate";
        case Benchmark::Mode::Increment:
            return "increment";
        case Benchmark::Mode::Push:
            return "push";
        case Benchmark::Mode::Pop:
            return "pop";
    }
}

[[nodiscard]] auto key_dist_name(Benchmark::KeyDistribution dist) {
    switch (dist) {
        case Benchmark::KeyDistribution::Random:
            return "uniform";
        case Benchmark::KeyDistribution::Ascending:
            return "ascending";
        case Benchmark::KeyDistribution::Descending:
            return "descending";
    }
    return "unknown";
}

[[nodiscard]] auto execution_mode_name(Benchmark::ExecutionMode exec_mode) {
    switch (exec_mode) {
        case Benchmark::ExecutionMode::Single:
            return "single";
        case Benchmark::ExecutionMode::Chunk:
            return "chunk";
        case Benchmark::ExecutionMode::Split:
            return "split";
    }
    return "unknown";
}

void write_settings(Benchmark::Settings const& settings, std::ostream& out) {
    out << "Threads: " << settings.num_threads << '\n'
        << "Prefill per thread: " << settings.prefill_per_thread << '\n'
        << "Operations per thread: " << settings.operations_per_thread << '\n'
        << "Mode: " << mode_name(settings.mode) << '\n';
    if (settings.mode == Benchmark::Mode::Increment) {
        out << "Increment range: [" << static_cast<unsigned long>(settings.min_increment) << ", "
            << static_cast<unsigned long>(settings.max_increment) << "]\n";
    } else if (settings.mode != Benchmark::Mode::Pop) {
        out << "Key distribution: " << key_dist_name(settings.key_distribution) << '\n';
        if (settings.key_distribution == Benchmark::KeyDistribution::Random) {
            out << "Key range: [" << static_cast<unsigned long>(settings.min_random) << ", "
                << static_cast<unsigned long>(settings.max_random) << "]\n";
        }
    }
    out << "Prefill range: [" << static_cast<unsigned long>(settings.min_prefill) << ", "
        << static_cast<unsigned long>(settings.max_prefill) << "]\n"
        << "Execution mode: " << execution_mode_name(settings.execution_mode) << '\n'
        << "Chunk size: " << settings.chunk_size << '\n'
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
        out << "Log operations to: " << settings_.log_file << '\n';
    }
#endif
#ifdef USE_UINT8
    out << "Data type: uint8_t\n";
#else
    out << "Data type: unsigned long\n";
#endif
}

int main(int argc, char* argv[]) {
    write_build_info(std::clog);
    Benchmark::Settings settings;
#ifdef QUALITY
    std::filesystem::path stats_file;
#endif
    cxxopts::Options cmd(argv[0]);
    cmd.add_options()
        // clang-format off
        ("h,help", "Print this help", cxxopts::value<bool>())
        ("j,threads", "Number of threads", cxxopts::value<int>(settings.num_threads), "NUMBER")
        ("p,prefill", "Prefill per thread", cxxopts::value<long long>(settings.prefill_per_thread), "NUMBER")
        ("n,ops", "Number of operations per thread", cxxopts::value<long long>(settings.operations_per_thread), "NUMBER")
        ("mode", "Operation mode ([a]lternate, [i]ncrement, p[u]sh, p[o]p)", cxxopts::value<char>(), "STRING")
        ("execution", "Execution mode ([c]hunk, [s]ingle, s[p]lit)", cxxopts::value<char>(), "STRING")
        ("key-dist", "Key distribution ([r]andom, [a]scending, [d]escending)", cxxopts::value<char>(), "STRING")
        ("min-prefill", "Min prefill key", cxxopts::value<key_type>(settings.min_prefill), "NUMBER")
        ("max-prefill", "Max prefill key", cxxopts::value<key_type>(settings.max_prefill), "NUMBER")
        ("min-increment", "Min increment", cxxopts::value<key_type>(settings.min_increment), "NUMBER")
        ("max-increment", "Max increment", cxxopts::value<key_type>(settings.max_increment), "NUMBER")
        ("min-random", "Min random key", cxxopts::value<key_type>(settings.min_random), "NUMBER")
        ("max-random", "Max random key", cxxopts::value<key_type>(settings.max_random), "NUMBER")
        ("chunk-size", "Chunk size", cxxopts::value<long long>(settings.chunk_size), "NUMBER")
        ("s,seed", "Initial seed", cxxopts::value<int>(settings.seed), "NUMBER")
#ifdef QUALITY
        ("stats-file", "File to write the stats to", cxxopts::value<std::filesystem::path>(stats_file), "PATH")
        ("l,log-file", "File to write the operation log to", cxxopts::value<std::filesystem::path>(settings.log_file), "PATH")
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
        if (args.count("mode") > 0) {
            settings.mode = parse_mode(args["mode"].as<char>());
        }
        if (args.count("key-dist") > 0) {
            settings.key_distribution = parse_key_distribution(args["key-dist"].as<char>());
        }
        if (args.count("execution") > 0) {
            settings.execution_mode = parse_execution_mode(args["execution"].as<char>());
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
    std::ofstream stats_out;
    if (!stats_file.empty()) {
        stats_out = std::ofstream(stats_file);
        if (!stats_out) {
            std::cerr << "Error: Could not open stats file " << stats_file << " for writing" << std::endl;
            return EXIT_FAILURE;
        }
    }
    std::ofstream log_out;
    if (!settings_.log_file.empty()) {
        log_out = std::ofstream(settings_.log_file);
        if (!log_out) {
            std::cerr << "Error: Could not open log file " << settings_.log_file << " for writing" << std::endl;
            return EXIT_FAILURE;
        }
    }
#endif

    Benchmark benchmark(settings);
    auto result = benchmark.run();
}
