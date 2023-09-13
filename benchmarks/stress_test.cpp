#include "build_info.hpp"
#include "priority_queue_selection.hpp"
#include "task.hpp"

#include "cxxopts.hpp"

#ifdef QUALITY
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
static constexpr auto max_key = static_cast<key_type>(1) << 63;
#endif

using pq_type = PriorityQueue<key_type, value_type, true>;
using handle_type = pq_type::handle_type;

enum class Mode { Increment, Alternate, Push, Pop };
enum class KeyDistribution { Random, Ascending, Descending };
enum class ExecutionMode { Chunk, Split };

[[nodiscard]] auto mode_name(Mode mode) {
    switch (mode) {
        case Mode::Alternate:
            return "alternate";
        case Mode::Increment:
            return "increment";
        case Mode::Push:
            return "push";
        case Mode::Pop:
            return "pop";
        default:
            return "";
    }
}

[[nodiscard]] auto key_dist_name(KeyDistribution dist) {
    switch (dist) {
        case KeyDistribution::Random:
            return "uniform";
        case KeyDistribution::Ascending:
            return "ascending";
        case KeyDistribution::Descending:
            return "descending";
        default:
            return "";
    }
}

[[nodiscard]] auto execution_mode_name(ExecutionMode exec_mode) {
    switch (exec_mode) {
        case ExecutionMode::Chunk:
            return "chunk";
        case ExecutionMode::Split:
            return "split";
        default:
            return "";
    }
}

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
    std::chrono::milliseconds timeout{0};
#ifdef QUALITY
    std::filesystem::path log_file;
    std::filesystem::path stats_file;
    bool per_element = false;
#endif
#ifdef USE_PAPI
    std::vector<std::string> papi_events;
#endif

    static void write(Settings const& settings, std::ostream& out) {
        out << "Threads: " << settings.num_threads << '\n'
            << "Prefill per thread: " << settings.prefill_per_thread << '\n'
            << "Operations per thread: " << settings.operations_per_thread << '\n'
            << "Mode: " << mode_name(settings.mode) << '\n';
        if (settings.mode == Mode::Increment) {
            out << "Increment range: [" << static_cast<unsigned long>(settings.min_increment) << ", "
                << static_cast<unsigned long>(settings.max_increment) << "]\n";
        } else if (settings.mode != Mode::Pop) {
            out << "Key distribution: " << key_dist_name(settings.key_distribution) << '\n';
            if (settings.key_distribution == KeyDistribution::Random) {
                out << "Key range: [" << static_cast<unsigned long>(settings.min_random) << ", "
                    << static_cast<unsigned long>(settings.max_random) << "]\n";
            }
        }
        out << "Prefill range: [" << static_cast<unsigned long>(settings.min_prefill) << ", "
            << static_cast<unsigned long>(settings.max_prefill) << "]\n"
            << "Execution mode: " << execution_mode_name(settings.execution_mode) << '\n';
        if (settings.execution_mode == ExecutionMode::Chunk) {
            out << "Chunk size: " << settings.chunk_size << '\n';
        }
        out << "Seed: " << settings.seed << '\n';
        if (settings.timeout.count() > 0) {
            out << "Timeout (ms): " << settings.timeout.count() << "ms\n";
        }
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
#endif
#ifdef USE_UINT8
        out << "Data type: uint8_t\n";
#else
        out << "Data type: unsigned long\n";
#endif
    }
};

struct Stats {
    std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::time_point::max();
    std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::time_point::min();
    long long num_pops = 0;
    long long num_pushes = 0;
    long long num_failed_pops = 0;
#ifdef USE_PAPI
    bool papi_success = true;
    std::vector<long long> papi_counter{};
#endif
#ifdef MQ_COUNT_STATS
    multiqueue::Counters mq_stats{};
#endif
    template <typename It>
    static Stats accumulate(It begin, It end) {
        Stats total_stats;
#ifdef USE_PAPI
        total_stats.papi_counter.resize(begin->papi_counter.size());
#endif
        for (; begin != end; ++begin) {
            total_stats.start_time = std::min(total_stats.start_time, begin->start_time);
            total_stats.end_time = std::max(total_stats.end_time, begin->end_time);
#ifdef USE_PAPI
            if (total_stats.papi_success) {
                if (begin->papi_success) {
                    std::transform(total_stats.papi_counter.begin(), total_stats.papi_counter.end(),
                                   begin->papi_counter.begin(), total_stats.papi_counter.begin(), std::plus<>());
                } else {
                    total_stats.papi_success = false;
                    std::fill(total_stats.papi_counter.begin(), total_stats.papi_counter.end(), 0);
                }
            }
#endif
            total_stats.num_pushes += begin->num_pushes;
            total_stats.num_pops += begin->num_pops;
            total_stats.num_failed_pops += begin->num_failed_pops;
#ifdef MQ_COUNT_STATS
            total_stats.mq_stats.locked_push_pq += begin->mq_stats.locked_push_pq;
            total_stats.mq_stats.empty_pop_pqs += begin->mq_stats.empty_pop_pqs;
            total_stats.mq_stats.locked_pop_pq += begin->mq_stats.locked_pop_pq;
            total_stats.mq_stats.stale_pop_pq += begin->mq_stats.stale_pop_pq;
#endif
        }
        return total_stats;
    }

    static void write(Stats const& stats, bool header, bool accumulated, std::ostream& out) {
        if (header) {
            out << "stat,time";
            out << ",pushes,pops,failed-pops";
#ifdef MQ_COUNT_STATS
            out << ",locked-push-pq,empty-pop-pqs,locked-pop-pq,stale-pop-pq";
#endif
#ifdef USE_PAPI
            for (std::size_t i = 0; i < stats.papi_counter.size(); ++i) {
                out << ",pc-" << i;
            }
#endif
            out << '\n';
        }
        out << (accumulated ? "total" : "single");
        out << ',' << (stats.end_time - stats.start_time).count();
        // clang-format off
        out << ',' << stats.num_pushes
            << ',' << stats.num_pops
            << ',' << stats.num_failed_pops;
        // clang-format on
#ifdef MQ_COUNT_STATS
        // clang-format off
        out << ',' << stats.mq_stats.locked_push_pq
            << ',' << stats.mq_stats.empty_pop_pqs
            << ',' << stats.mq_stats.locked_pop_pq
            << ',' << stats.mq_stats.stale_pop_pq;
        // clang-format on
#endif
#ifdef USE_PAPI
        for (auto c : stats.papi_counter) {
            out << ',' << c;
        }
#endif
        out << '\n';
    }
};

struct ThreadData {
    task::Control tc;
    handle_type handle;
#ifdef QUALITY
    operation_log::OperationLog op_log{};
#endif
    Stats stats{};
};

int pushes_per_iteration(Mode mode) {
    switch (mode) {
        case Mode::Alternate:
        case Mode::Increment:
        case Mode::Push:
            return 1;
        default:
            return 0;
    }
}

int pops_per_iteration(Mode mode) {
    switch (mode) {
        case Mode::Alternate:
        case Mode::Increment:
        case Mode::Pop:
            return 1;
        default:
            return 0;
    }
}

long long num_iterations(Settings const& settings) {
    return settings.num_threads * settings.operations_per_thread /
        (pushes_per_iteration(settings.mode) + pops_per_iteration(settings.mode));
}

long long num_pushes(Settings const& settings) {
    return settings.num_threads * settings.prefill_per_thread +
        pushes_per_iteration(settings.mode) * num_iterations(settings);
}

long long num_pops(Settings const& settings) {
    return pops_per_iteration(settings.mode) * num_iterations(settings);
}

long long max_num_elements(Settings const& settings) {
    return settings.num_threads *
        (settings.prefill_per_thread + (settings.mode == Mode::Push ? settings.operations_per_thread : 0));
}

bool validate(Settings const& settings) {
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
        std::cerr << "Error: Max must be greater than min\n";
        return false;
    }
    if (settings.chunk_size <= 0) {
        std::cerr << "Error: Chunk size must be greater than 0\n";
        return false;
    }
    if (settings.prefill_per_thread < 0 || settings.operations_per_thread < 0) {
        std::cerr << "Error: Prefill and operations must be nonnegative\n";
        return false;
    }
    if (num_pops(settings) > num_pushes(settings)) {
        std::cerr << "Error: Number of pushes must not be less than number of pops\n";
        return false;
    }
#ifdef QUALITY
    if (num_pops(settings) == 0) {
        std::cerr << "Error: Number of pops must be greater than 0 to measure quality\n";
        return false;
    }
    assert(num_pushes(settings) > 0);
    if (static_cast<value_type>(num_pushes(settings) - 1) > std::numeric_limits<value_type>::max()) {
        std::cerr << "Error: Number of pushes must be representable in `value_type`\n";
        return false;
    }
#endif
#ifdef USE_PAPI
    for (auto const& name : settings.papi_events) {
        if (PAPI_query_named_event(name.c_str()) != PAPI_OK) {
            std::cerr << "Error: PAPI event " << name << " not available\n";
            return false;
        }
    }
#endif
    return true;
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

void generate_keys(Settings const& settings, std::vector<key_type>& keys, int id) {
    std::seed_seq seed{settings.seed, id};
    std::default_random_engine rng(seed);

    auto prefill_begin = keys.begin() + id * settings.prefill_per_thread;
    auto prefill_end = prefill_begin + settings.prefill_per_thread;
    auto prefill_dist = std::uniform_int_distribution<key_type>(settings.min_prefill, settings.max_prefill);
    std::generate(prefill_begin, prefill_end, [&prefill_dist, &rng]() { return prefill_dist(rng); });
    auto first_workload = keys.begin() + settings.num_threads * settings.prefill_per_thread;
    auto workload_begin = first_workload + id * (std::distance(first_workload, keys.end()) / settings.num_threads);
    auto workload_end = first_workload + (id + 1) * (std::distance(first_workload, keys.end()) / settings.num_threads);
    switch (settings.mode) {
        case Mode::Increment: {
            auto dist = std::uniform_int_distribution<key_type>(settings.min_increment, settings.max_increment);
            std::generate(workload_begin, workload_end, [&dist, &rng]() { return dist(rng); });
            break;
        }
        case Mode::Push:
        case Mode::Alternate: {
            switch (settings.key_distribution) {
                case KeyDistribution::Random: {
                    auto dist = std::uniform_int_distribution<key_type>(settings.min_random, settings.max_random);
                    std::generate(workload_begin, workload_end, [&dist, &rng]() { return dist(rng); });
                    break;
                }
                case KeyDistribution::Ascending: {
                    std::iota(workload_begin, workload_end, std::distance(keys.begin(), workload_begin));
                    break;
                }
                case KeyDistribution::Descending: {
                    std::iota(std::reverse_iterator(workload_end), std::reverse_iterator(workload_begin),
                              std::distance(workload_end, keys.end()));
                    break;
                }
            }
        }
        case Mode::Pop:
            break;
    }
}

#ifdef QUALITY
void push(key_type key, std::size_t index, ThreadData& data) {
    data.handle.push({key, static_cast<value_type>(index)});
    auto tick = operation_log::get_tick();
    data.op_log.pushes.push_back({tick, key, index});
    ++data.stats.num_pushes;
}

key_type pop(ThreadData& data) {
    do {
        auto tick = operation_log::get_tick();
        auto retval = data.handle.try_pop();
        if (retval) {
            data.op_log.pops.push_back({tick, static_cast<std::size_t>(retval->second)});
            ++data.stats.num_pops;
            return retval->first;
        }
        ++data.stats.num_failed_pops;
    } while (true);
}
#else
void push(key_type key, std::size_t index, ThreadData& data) {
    data.handle.push({key, static_cast<value_type>(index)});
    ++data.stats.num_pushes;
};

auto pop(ThreadData& data) {
    do {
        auto retval = data.handle.try_pop();
        if (retval) {
            ++data.stats.num_pops;
            return retval->first;
        }
        ++data.stats.num_failed_pops;
    } while (true);
};
#endif

void prefill(Settings const& settings, std::vector<key_type> const& keys, ThreadData& data) {
    auto offset = static_cast<std::size_t>(data.tc.id() * settings.prefill_per_thread);
    for (std::size_t i = 0; i < static_cast<std::size_t>(settings.prefill_per_thread); ++i) {
        push(keys[offset + i], offset + i, data);
    }
}

template <typename Executor>
void execute_benchmark(Settings const& settings, std::vector<key_type> const& keys, ThreadData& data,
                       Executor& executor) {
    auto offset = static_cast<std::size_t>(settings.num_threads * settings.prefill_per_thread);
#ifdef MQ_COUNT_STATS
    data.handle.reset_counters();
#endif
    switch (settings.mode) {
        case Mode::Alternate: {
            auto work = [&, offset](std::size_t index) {
                push(keys[offset + index], offset + index, data);
                pop(data);
            };
            execute(settings, executor, work, data);
            break;
        }
        case Mode::Increment: {
            auto work = [&, offset](std::size_t index) {
                auto key = pop(data);
                push(key + keys[offset + index], offset + index, data);
            };
            execute(settings, executor, work, data);
            break;
        }
        case Mode::Push: {
            auto work = [&, offset](std::size_t index) { push(keys[offset + index], offset + index, data); };
            execute(settings, executor, work, data);
            break;
        }
        case Mode::Pop: {
            auto work = [&](std::size_t) { pop(data); };
            execute(settings, executor, work, data);
            break;
        }
    }
#ifdef MQ_COUNT_STATS
    data.stats.mq_stats = data.handle.get_counters();
#endif
}

template <typename Executor, typename Work>
void execute(Settings const& settings, Executor& executor, Work work, ThreadData& data) {
#ifdef USE_PAPI
    int event_set = PAPI_NULL;
    bool papi_started = false;
    if (!settings.papi_events.empty()) {
        data.stats.papi_counter.resize(settings.papi_events.size());
        try {
            event_set = start_papi(settings.papi_events);
            papi_started = true;
        } catch (std::runtime_error const& e) {
            data.tc.write(std::cerr) << "Error: " << e.what() << '\n';
            data.stats.papi_success = false;
        }
    }
#endif
    auto [t_start, t_end] = executor(data.tc, static_cast<std::size_t>(num_iterations(settings)), work);

#ifdef USE_PAPI
    if (papi_started) {
        if (int ret = PAPI_stop(event_set, data.stats.papi_counter.data()); ret != PAPI_OK) {
            data.tc.write(std::cerr) << "Error: Failed to stop performance counters\n";
            data.stats.papi_success = false;
        }
    }
#endif
    data.stats.start_time = t_start;
    data.stats.end_time = t_end;
}

template <typename Executor>
void run_thread(Settings const& settings, std::vector<key_type>& keys, ThreadData& data, Executor& executor) {
#ifdef QUALITY
    data.op_log.pushes.reserve(static_cast<std::size_t>(num_pushes(settings)));
    data.op_log.pops.reserve(static_cast<std::size_t>(num_pops(settings)));
#endif
    std::chrono::steady_clock::time_point t_start;
    std::chrono::steady_clock::time_point t_end;
    data.tc.once([]() { std::clog << "Generating keys..." << std::flush; });
    data.tc.synchronize();
    data.tc.once([&t_start]() { t_start = std::chrono::steady_clock::now(); });
    data.tc.synchronize();
    generate_keys(settings, keys, data.tc.id());
    data.tc.synchronize();
    data.tc.once([&t_end, &t_start]() {
        t_end = std::chrono::steady_clock::now();
        std::clog << "done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count() << "ms)"
                  << std::endl;
        std::clog << "Prefilling..." << std::flush;
        t_start = std::chrono::steady_clock::now();
    });
    data.tc.synchronize();
    prefill(settings, keys, data);
    data.tc.synchronize();
    data.tc.once([&t_end, &t_start]() {
        t_end = std::chrono::steady_clock::now();
        std::clog << "done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count() << "ms)"
                  << std::endl;
        std::clog << "Running benchmark..." << std::flush;
        t_start = std::chrono::steady_clock::now();
    });
    execute_benchmark(settings, keys, data, executor);
    data.tc.once([&t_end, &t_start]() {
        t_end = std::chrono::steady_clock::now();
        std::clog << "done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count() << "ms)"
                  << std::endl;
    });
}

template <typename Integer = std::size_t>
class ChunkExecutor {
    std::atomic<Integer> counter_;
    Integer chunk_size_;
    std::chrono::milliseconds timeout_;

   public:
    explicit ChunkExecutor(Integer chunk_size, std::chrono::milliseconds t = std::chrono::milliseconds{0})
        : chunk_size_{chunk_size}, timeout_{t} {
    }

    template <typename Work, typename... Args>
    auto operator()(task::Control const& tc, Integer max, Work work, Args&&... args) {
        if (tc.id() == 0) {
            counter_.store(0, std::memory_order_relaxed);
        }
        tc.synchronize();
        auto t_start = std::chrono::high_resolution_clock::now();
        for (Integer from = counter_.fetch_add(chunk_size_, std::memory_order_relaxed); from < max;
             from = counter_.fetch_add(chunk_size_, std::memory_order_relaxed)) {
            if (timeout_.count() > 0) {
                auto t_now = std::chrono::high_resolution_clock::now();
                if (t_now - t_start > timeout_) {
                    break;
                }
            }
            auto to = std::min(from + chunk_size_, max);
            for (auto i = from; i < to; ++i) {
                work(i, args...);  // no perfect forwarding here
            }
        }
        auto t_end = std::chrono::high_resolution_clock::now();
        tc.synchronize();
        return std::pair{t_start, t_end};
    }
};

template <typename Integer = std::size_t>
class SplitExecutor {
    std::chrono::milliseconds timeout_;

   public:
    explicit SplitExecutor(std::chrono::milliseconds t = std::chrono::milliseconds{0}) : timeout_{t} {
    }
    template <typename Work, typename... Args>
    auto operator()(task::Control const& tc, Integer max, Work work, Args&&... args) {
        auto from = (static_cast<Integer>(tc.id()) * max) / static_cast<Integer>(tc.num_threads());
        auto to = tc.id() == tc.num_threads() - 1 ? max : from + (max / static_cast<Integer>(tc.num_threads()));
        tc.synchronize();
        auto t_start = std::chrono::high_resolution_clock::now();
        for (; from != to; ++from) {
            if (timeout_.count() > 0) {
                auto t_now = std::chrono::high_resolution_clock::now();
                if (t_now - t_start > timeout_) {
                    break;
                }
            }
            work(from, args...);  // no perfect forwarding here
        }
        auto t_end = std::chrono::high_resolution_clock::now();
        tc.synchronize();
        return std::pair{t_start, t_end};
    }
};

static void write_all_stats(std::vector<Stats> const& stats, std::ostream& out) {
    Stats::write(Stats::accumulate(stats.begin(), stats.end()), true, true, out);
    for (auto const& s : stats) {
        Stats::write(s, false, false, out);
    }
}

bool run_benchmark(Settings const& settings, pq_type& pq) {
#ifdef QUALITY
    std::ofstream stats_out;
    if (!settings.stats_file.empty()) {
        stats_out = std::ofstream(settings.stats_file);
        if (!stats_out) {
            std::cerr << "Error: Could not open stats file " << settings.stats_file << " for writing" << std::endl;
            return false;
        }
    }
    std::ofstream log_out;
    if (!settings.log_file.empty()) {
        log_out = std::ofstream(settings.log_file);
        if (!log_out) {
            std::cerr << "Error: Could not open log file " << settings.log_file << " for writing" << std::endl;
            return false;
        }
    }
    auto op_logs = std::vector<operation_log::OperationLog>(static_cast<std::size_t>(settings.num_threads));
#endif

    auto keys = std::vector<key_type>(static_cast<std::size_t>(num_pushes(settings)));
    auto stats = std::vector<Stats>(static_cast<std::size_t>(settings.num_threads));

    auto run_with_executor = [&](auto& executor) {
        affinity::NUMA numa_affinity{cores_per_numa_node, num_numa_nodes};
        task::Runner runner(numa_affinity, settings.num_threads, [&](auto tc) {
            auto data = ThreadData{tc, pq.get_handle()};
            run_thread(settings, keys, data, executor);
            stats[static_cast<std::size_t>(tc.id())] = std::move(data.stats);
#ifdef QUALITY
            op_logs[static_cast<std::size_t>(tc.id())] = std::move(data.op_log);
#endif
        });
        runner.wait();
    };

    switch (settings.execution_mode) {
        case ExecutionMode::Chunk: {
            ChunkExecutor executor{static_cast<std::size_t>(settings.chunk_size), settings.timeout};
            run_with_executor(executor);
            break;
        }
        case ExecutionMode::Split: {
            SplitExecutor executor{settings.timeout};
            run_with_executor(executor);
            break;
        }
    }

#ifdef QUALITY
    if (!settings.stats_file.empty()) {
        write_all_stats(stats, stats_out);
        stats_out.close();
    }
    std::clog << "Merging logs..." << std::flush;
    auto t_start = std::chrono::steady_clock::now();
    auto merged_logs = operation_log::merge_logs(op_logs);
    auto t_end = std::chrono::steady_clock::now();
    std::clog << "done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count() << "ms)\n";
    if (!settings.log_file.empty()) {
        std::clog << "Writing operation logs..." << std::flush;
        t_start = std::chrono::steady_clock::now();
        operation_log::write_logs(merged_logs, log_out);
        log_out.close();
        t_end = std::chrono::steady_clock::now();
        std::clog << "done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count()
                  << "ms)\n";
    }
    if (!operation_log::verify_logs(merged_logs)) {
        return false;
    }
    std::clog << "Replaying logs..." << std::flush;
    t_start = std::chrono::steady_clock::now();
    auto metrics = operation_log::replay_logs(std::move(merged_logs));
    t_end = std::chrono::steady_clock::now();
    std::clog << "done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count() << "ms)\n";
    if (settings.per_element) {
        operation_log::write_metrics(metrics, std::cout);
    } else {
        operation_log::write_metrics_average(metrics, std::cout);
    }
#else
    write_all_stats(stats, std::cout);
#endif
    return true;
}

Mode parse_mode(char c) {
    switch (c) {
        case 'a':
            return Mode::Alternate;
        case 'i':
            return Mode::Increment;
        case 'o':
            return Mode::Pop;
        case 'u':
            return Mode::Push;
        default:
            throw std::invalid_argument("Invalid work mode");
    }
}

KeyDistribution parse_key_distribution(char c) {
    switch (c) {
        case 'r':
            return KeyDistribution::Random;
        case 'a':
            return KeyDistribution::Ascending;
        case 'd':
            return KeyDistribution::Descending;
        default:
            throw std::invalid_argument("Invalid key distribution");
    }
}

ExecutionMode parse_execution_mode(char c) {
    switch (c) {
        case 'c':
            return ExecutionMode::Chunk;
        case 's':
            return ExecutionMode::Split;
        default:
            throw std::invalid_argument("Invalid execution mode");
    }
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
        ("timeout", "Timeout in milliseconds", cxxopts::value<int>(timeout_ms), "NUMBER")
#ifdef QUALITY
        ("stats-file", "File to write the stats to", cxxopts::value<std::filesystem::path>(settings.stats_file), "PATH")
        ("l,log-file", "File to write the operation log to", cxxopts::value<std::filesystem::path>(settings.log_file), "PATH")
        ("per-element", "Output metrics per element", cxxopts::value<bool>(settings.per_element), "NUMBER")
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
        settings.timeout = std::chrono::milliseconds(timeout_ms);
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
    Settings::write(settings, std::clog);

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

    if (!validate(settings)) {
        return EXIT_FAILURE;
    }

    auto pq = pq_type(settings.num_threads, static_cast<std::size_t>(num_pushes(settings)), settings.pq_settings);
    std::clog << "Priority queue: ";
    pq.describe(std::clog) << '\n';

    bool success = run_benchmark(settings, pq);
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
