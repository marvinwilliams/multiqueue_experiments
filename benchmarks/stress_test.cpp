#include "build_info.hpp"
#include "task.hpp"
#include "wrapper/priority.hpp"
#include "wrapper/selector.hpp"

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

using pq_type = PriorityQueue<key_type, value_type, Priority::Min>;
using handle_type = pq_type::handle_type;

struct Settings {
    enum class Mode { Update, Random, Pop, PushRandom, PushAscending };

    int num_threads = 4;
    long long prefill_per_thread = 1 << 20;
    long long iterations_per_thread = 1 << 24;
    Mode mode = Mode::Update;
    pq_type::config_type pq_settings{};
    key_type min_prefill = 1;
    key_type max_prefill = 1 << 20;
    long min_update = 0;
    long max_update = 1 << 20;
    key_type min_random = 1;
    key_type max_random = 1 << 20;
    long long chunk_size = 1 << 12;
    int seed = 1;
    std::chrono::seconds timeout{0};
#ifdef QUALITY
    std::filesystem::path log_file;
    std::filesystem::path metric_file;
#endif
    std::filesystem::path thread_stat_file;
#ifdef USE_PAPI
    std::vector<std::string> papi_events;
#endif

    static constexpr int pushes_per_iteration(Mode mode) noexcept {
        return mode == Mode::Pop ? 0 : 1;
    }

    static constexpr int pops_per_iteration(Mode mode) noexcept {
        return mode == Mode::PushRandom || mode == Mode::PushAscending ? 0 : 1;
    }

    static void write(Settings const& settings, std::ostream& out) {
        auto mode_name = [mode = settings.mode]() {
            switch (mode) {
                case Mode::Update:
                    return "update";
                case Mode::Random:
                    return "random";
                case Mode::Pop:
                    return "pop";
                case Mode::PushRandom:
                    return "push (random)";
                case Mode::PushAscending:
                    return "push (ascending)";
            }
        };

        out << "Threads: " << settings.num_threads << '\n'
            << "Prefill per thread: " << settings.prefill_per_thread << '\n'
            << "Prefill range: [" << static_cast<unsigned long>(settings.min_prefill) << ", "
            << static_cast<unsigned long>(settings.max_prefill) << "]\n"
            << "Iterations per thread: " << settings.iterations_per_thread << '\n'
            << "Mode: " << mode_name() << '\n';
        if (settings.mode == Mode::Update) {
            out << "Update range: [" << settings.min_update << ", " << settings.max_update << "]\n";
        } else if (settings.mode == Mode::PushRandom || settings.mode == Mode::Random) {
            out << "Key range: [" << static_cast<unsigned long>(settings.min_random) << ", "
                << static_cast<unsigned long>(settings.max_random) << "]\n";
        }
        out << "Chunk size: " << settings.chunk_size << '\n' << "Seed: " << settings.seed << '\n';
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
        if (!settings.metric_file.empty()) {
            out << "Write per element metrics to: " << settings.metric_file << '\n';
        }
#endif
        if (!settings.thread_stat_file.empty()) {
            out << "Write thread stats to: " << settings.thread_stat_file << '\n';
        }
#ifdef USE_UINT8
        out << "Data type: uint8_t\n";
#else
        out << "Data type: unsigned long\n";
#endif
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
            settings.max_update < settings.min_update) {
            std::cerr << "Error: Max must be greater than min\n";
            return false;
        }
        if (settings.chunk_size <= 0) {
            std::cerr << "Error: Chunk size must be greater than 0\n";
            return false;
        }
        if (settings.prefill_per_thread < 0 || settings.iterations_per_thread < 0) {
            std::cerr << "Error: Prefill and iterations must be nonnegative\n";
            return false;
        }
        if (settings.iterations_per_thread * pops_per_iteration(settings.mode) >
            settings.prefill_per_thread + settings.iterations_per_thread * pushes_per_iteration(settings.mode)) {
            std::cerr << "Error: Number of pops must not be greater than number of pushes\n";
            return false;
        }
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
};

struct BenchmarkData {
    std::vector<key_type> prefill;
    std::vector<long long> push_data;
};

void generate_benchmark_data(Settings const& settings, task::Control const& tc, BenchmarkData& benchmark_data) {
    std::seed_seq seed{settings.seed, tc.id()};
    std::default_random_engine rng(seed);

    auto prefill_begin = benchmark_data.prefill.begin() + tc.id() * settings.prefill_per_thread;
    std::generate(prefill_begin, prefill_begin + settings.prefill_per_thread,
                  [&rng, min = settings.min_prefill, max = settings.max_prefill]() {
                      return std::uniform_int_distribution<key_type>(min, max)(rng);
                  });
    auto push_data_begin = benchmark_data.push_data.begin() +
        tc.id() * settings.iterations_per_thread * Settings::pushes_per_iteration(settings.mode);
    switch (settings.mode) {
        case Settings::Mode::Update: {
            std::generate(
                push_data_begin,
                push_data_begin + settings.iterations_per_thread * Settings::pushes_per_iteration(settings.mode),
                [&rng, min = settings.min_update, max = settings.max_update]() {
                    return std::uniform_int_distribution<long>(min, max)(rng);
                });
            break;
        }
        case Settings::Mode::Random:
        case Settings::Mode::PushRandom: {
            auto dist = std::uniform_int_distribution<key_type>(settings.min_random, settings.max_random);
            std::generate(
                push_data_begin,
                push_data_begin + settings.iterations_per_thread * Settings::pushes_per_iteration(settings.mode),
                [&rng, min = settings.min_random, max = settings.max_random]() {
                    return static_cast<long>(std::uniform_int_distribution<key_type>(min, max)(rng));
                });
            break;
        }
        case Settings::Mode::PushAscending: {
            std::iota(push_data_begin,
                      push_data_begin + settings.iterations_per_thread * Settings::pushes_per_iteration(settings.mode),
                      std::distance(benchmark_data.push_data.begin(), push_data_begin));
            break;
        }
        case Settings::Mode::Pop:
            break;
    }
}

class Executor {
    std::atomic<std::size_t> counter_;

   public:
    struct Result {
        std::chrono::high_resolution_clock::duration work_time;
        long long work_count;
    };

    template <typename Work, typename... Args>
    Result operator()(task::Control const& tc, std::size_t max, std::size_t chunk_size, std::chrono::seconds timeout,
                      Work work, Args&&... args) {
        if (tc.id() == 0) {
            counter_.store(0, std::memory_order_relaxed);
        }
        Result result{};
        tc.synchronize();
        auto start_time = std::chrono::high_resolution_clock::now();
        for (auto from = counter_.fetch_add(chunk_size, std::memory_order_relaxed); from < max;
             from = counter_.fetch_add(chunk_size, std::memory_order_relaxed)) {
            auto to = std::min(from + chunk_size, max);
            for (auto i = from; i < to; ++i) {
                work(i, args...);  // no perfect forwarding here
            }
            result.work_count += static_cast<long long>(to - from);
            if (timeout.count() > 0) {
                auto now = std::chrono::high_resolution_clock::now();
                if (now - start_time > timeout) {
                    break;
                }
            }
        }
        result.work_time = std::chrono::high_resolution_clock::now() - start_time;
        tc.synchronize();
        return result;
    }
};

struct Context {
    pq_type pq;
    BenchmarkData benchmark_data;
    Settings settings;
    Executor executor;
};

struct Stats {
    std::chrono::high_resolution_clock::duration work_time{};
    long long num_iterations = 0;
    long long num_failed_pops = 0;
#ifdef USE_PAPI
    bool papi_success = true;
    std::vector<long long> papi_counter{};
#endif
#ifdef MQ_COUNT_STATS
    multiqueue::Counters mq_stats{};
#endif
#ifdef QUALITY
    operation_log::OperationLog op_log{};
#endif

    template <typename It>
    static Stats accumulate(It begin, It end) {
        if (begin == end) {
            return {};
        }
        Stats accum_stats = *begin;
#ifdef QUALITY
        accum_stats.op_log.pushes.reserve(
            std::accumulate(begin, end, std::size_t{0},
                            [](std::size_t sum, auto const& stat) { return sum + stat.op_log.pushes.size(); }));
        accum_stats.op_log.pops.reserve(
            std::accumulate(logs.begin(), logs.end(), std::size_t{0},
                            [](std::size_t sum, auto const& stat) { return sum + stat.op_log.pops.size(); }));
#endif
        ++begin;
        for (; begin != end; ++begin) {
            accum_stats.work_time += begin->work_time;
#ifdef USE_PAPI
            accum_stats.papi_success &= begin->papi_success;
            if (accum_stats.papi_success) {
                std::transform(accum_stats.papi_counter.begin(), accum_stats.papi_counter.end(),
                               begin->papi_counter.begin(), accum_stats.papi_counter.begin(), std::plus<>());
            }
#endif
            accum_stats.num_iterations += begin->num_iterations;
            accum_stats.num_failed_pops += begin->num_failed_pops;
#ifdef MQ_COUNT_STATS
            accum_stats.mq_stats.locked_push_pq += begin->mq_stats.locked_push_pq;
            accum_stats.mq_stats.empty_pop_pqs += begin->mq_stats.empty_pop_pqs;
            accum_stats.mq_stats.locked_pop_pq += begin->mq_stats.locked_pop_pq;
            accum_stats.mq_stats.stale_pop_pq += begin->mq_stats.stale_pop_pq;
#endif
#ifdef QUALITY
            accum_stats.op_log.pushes.insert(accum_stats.op_log.pushes.end(), begin->op_log.pushes.begin(),
                                             begin->op_log.pushes.end());
            accum_stats.op_log.pops.insert(accum_stats.op_log.pops.end(), begin->op_log.pops.begin(),
                                           begin->op_log.pops.end());
#endif
        }
#ifdef QUALITY
        std::sort(accum_stats.op_log.pushes.begin(), accum_stats.op_log.pushes.end());
        std::sort(accum_stats.op_log.pops.begin(), accum_stats.op_log.pops.end());
#endif
        return accum_stats;
    }

    static void write_header(Settings const& settings, std::ostream& out) {
        out << "time,iterations,failed-pops";
#ifdef MQ_COUNT_STATS
        out << ",locked-push-pq,empty-pop-pqs,locked-pop-pq,stale-pop-pq";
#endif
#ifdef USE_PAPI
        for (auto const& name : settings.papi_events) {
            out << "," << name;
        }
#endif
        out << '\n';
    }

    static void write(Stats const& stats, std::ostream& out) {
        out << stats.work_time.count() << ',' << stats.num_iterations << ',' << stats.num_iterations << ','
            << stats.num_failed_pops;
#ifdef MQ_COUNT_STATS
        out << ',' << stats.mq_stats.locked_push_pq << ',' << stats.mq_stats.empty_pop_pqs << ','
            << stats.mq_stats.locked_pop_pq << ',' << stats.mq_stats.stale_pop_pq;
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
    Stats stats{};
};

void pq_push(key_type key, value_type value, ThreadData& data) {
    data.handle.push({key, value});
#ifdef QUALITY
    auto tick = static_cast<long long>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    data.op_log.pushes.push_back({tick, key, static_cast<std::size_t>(value)});
#endif
}

key_type pq_pop(ThreadData& data) {
    while (true) {
#ifdef QUALITY
        auto tick = static_cast<long long>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
#endif
        auto retval = data.handle.try_pop();
        if (retval) {
#ifdef QUALITY
            data.op_log.pops.push_back({tick, static_cast<std::size_t>(retval->second)});
#endif
            return retval->first;
        }
        ++data.stats.num_failed_pops;
    }
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

void prefill(Settings const& settings, std::vector<key_type> const& prefill, ThreadData& data) {
    auto offset = static_cast<std::size_t>(data.tc.id() * settings.prefill_per_thread);
    for (auto i = offset; i < offset + static_cast<std::size_t>(settings.prefill_per_thread); ++i) {
        pq_push(prefill[i], static_cast<value_type>(i), data);
    }
}

void dispatch_benchmark_mode(Settings const& settings, std::vector<long> const& push_data, ThreadData& data) {
#ifdef MQ_COUNT_STATS
    data.handle.reset_counters();
#endif
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
    auto offset = static_cast<std::size_t>(settings.num_threads * settings.prefill_per_thread);
    switch (settings.mode) {
        case Settings::Mode::Update: {
            auto work = [&push_data, &data, offset](std::size_t i) {
                auto key = pq_pop(data);
                auto new_key = static_cast<key_type>(static_cast<long long>(key) + push_data[i]);
                pq_push(new_key, static_cast<value_type>(offset + i), data);
            };
            auto result = executor(data.tc, static_cast<std::size_t>(Settings::num_iterations(settings)), work);
            break;
        }
        case Settings::Mode::Random: {
            auto work = [&push_data, &data, offset](std::size_t i) {
                pq_pop(data);
                auto new_key = static_cast<key_type>(push_data[i]);
                pq_push(new_key, static_cast<value_type>(offset + i), data);
            };
            benchmark_mode(settings, work, data);
            break;
        }
        case Settings::Mode::PushRandom: {
            case Settings::Mode::PushAscending: {
                auto work = [&push_data, &data, offset](std::size_t i) {
                    auto new_key = static_cast<key_type>(push_data[i]);
                    pq_push(new_key, static_cast<value_type>(offset + i), data);
                };
                benchmark_mode(settings, work, data);
                break;
            }
            case Settings::Mode::Pop: {
                auto work = [&data](std::size_t) { pq_pop(data); };
                benchmark_mode(settings, work, data);
                break;
            }
        }
    }
#ifdef USE_PAPI
    if (papi_started) {
        if (int ret = PAPI_stop(event_set, data.stats.papi_counter.data()); ret != PAPI_OK) {
            data.tc.write(std::cerr) << "Error: Failed to stop performance counters\n";
            data.stats.papi_success = false;
        }
    }
#endif
    data.stats.start_time = result.start_time;
    data.stats.end_time = result.end_time;
    data.stats.num_pushes = result.work_count * Settings::pushes_per_iteration(settings.mode);
    data.stats.num_pops = result.work_count * Settings::pops_per_iteration(settings.mode);
#ifdef MQ_COUNT_STATS
    data.stats.mq_stats = data.handle.get_counters();
#endif
}

template <typename Executor>
Stats benchmark_thread(Settings const& settings, BenchmarkData& workload, ThreadData& data, Executor& executor) {
#ifdef QUALITY
    data.op_log.pushes.reserve(
        static_cast<std::size_t>(settings.num_threads * settings.prefill_per_thread + Settings::num_pushes(settings)));
    data.op_log.pops.reserve(static_cast<std::size_t>(Settings::num_pops(settings)));
#endif
    std::chrono::steady_clock::time_point t_start;
    std::chrono::steady_clock::time_point t_end;
    data.tc.once([]() { std::clog << "Generating workload..." << std::flush; });
    data.tc.synchronize();
    data.tc.once([&t_start]() { t_start = std::chrono::steady_clock::now(); });
    data.tc.synchronize();
    generate_benchmark_data(settings, workload, data.tc.id());
    data.tc.synchronize();
    data.tc.once([&t_end, &t_start]() {
        t_end = std::chrono::steady_clock::now();
        std::clog << "done (" << std::fixed << std::setprecision(2)
                  << std::chrono::duration<double>(t_end - t_start).count() << "s)" << std::endl;
        std::clog << "Prefilling..." << std::flush;
        t_start = std::chrono::steady_clock::now();
    });
    data.tc.synchronize();
    prefill(settings, workload.prefill, data);
    data.tc.synchronize();
    data.tc.once([&t_end, &t_start]() {
        t_end = std::chrono::steady_clock::now();
        std::clog << "done (" << std::chrono::duration<double>(t_end - t_start).count() << "s)" << std::endl;
        std::clog << "Running benchmark..." << std::flush;
        t_start = std::chrono::steady_clock::now();
    });
    execute_benchmark(settings, workload.push_data, data, executor);
    data.tc.once([&t_end, &t_start]() {
        t_end = std::chrono::steady_clock::now();
        std::clog << "done (" << std::chrono::duration<double>(t_end - t_start).count() << "s)" << std::endl;
    });
    return data.stats;
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

    BenchmarkData workload;
    workload.prefill.resize(static_cast<std::size_t>(settings.num_threads * settings.prefill_per_thread));
    workload.push_data.resize(static_cast<std::size_t>(Settings::num_pushes(settings)));
    auto stats = std::vector<Stats>(static_cast<std::size_t>(settings.num_threads));
    auto executor = Executor{static_cast<std::size_t>(settings.chunk_size), settings.timeout};
    affinity::NUMA numa_affinity{cores_per_numa_node, num_numa_nodes};
    task::Runner runner(numa_affinity, settings.num_threads, [&](auto tc) {
        auto data = ThreadData{tc, pq.get_handle()};
        stats[static_cast<std::size_t>(tc.id())] = run_thread(settings, workload, data, executor);
#ifdef QUALITY
        op_logs[static_cast<std::size_t>(tc.id())] = std::move(data.op_log);
#endif
    });
    runner.wait();

#ifdef QUALITY
    if (!settings.stats_file.empty()) {
        Stats::write(stats.begin(), stats.end(), settings, stats_out);
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
        out << "rank_error,delay\n";
        for (auto const& m : metrics) {
            out << m.rank_error << ',' << m.delay << '\n';
        }
    } else {
        assert(!metrics.empty());
        auto summed_metrics = std::reduce(metrics.begin(), metrics.end(), Metrics{}, [](Metrics a, Metrics const& b) {
            a.rank_error += b.rank_error;
            a.delay += b.delay;
            return a;
        });
        out << "rank_error,delay\n";
        out << std::fixed << std::setprecision(2)
            << static_cast<double>(summed_metrics.rank_error) / static_cast<double>(metrics.size()) << ','
            << static_cast<double>(summed_metrics.delay) / static_cast<double>(metrics.size()) << '\n';
    }
#else
    Stats::write(stats.begin(), stats.end(), settings, std::cout);
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
        ("execution", "Execution mode ([c]hunk, [s]plit)", cxxopts::value<char>(), "STRING")
        ("key-dist", "Key distribution ([r]andom, [a]scending, [d]escending)", cxxopts::value<char>(), "STRING")
        ("min-prefill", "Min prefill key", cxxopts::value<key_type>(settings.min_prefill), "NUMBER")
        ("max-prefill", "Max prefill key", cxxopts::value<key_type>(settings.max_prefill), "NUMBER")
        ("min-increment", "Min increment", cxxopts::value<long>(settings.min_increment), "NUMBER")
        ("max-increment", "Max increment", cxxopts::value<long>(settings.max_increment), "NUMBER")
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
        ("per-thread", "Output metrics per thread", cxxopts::value<bool>(settings.per_thread), "NUMBER")
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

    if (!Settings::validate(settings)) {
        return EXIT_FAILURE;
    }

    auto max_num_elements = settings.num_threads *
        (settings.prefill_per_thread +
         (settings.mode == Mode::PushRandom || settings.mode == Mode::PushAscending ? settings.iterations_per_thread
                                                                                    : 0));

    auto pq = pq_type(settings.num_threads, static_cast<std::size_t>(2 * max_num_elements), settings.pq_settings);
    std::clog << "Priority queue: ";
    pq.describe(std::clog) << '\n';

    bool success = run_benchmark(settings, pq);
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
