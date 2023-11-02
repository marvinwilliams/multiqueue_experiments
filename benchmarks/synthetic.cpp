#include "build_info.hpp"
#include "task.hpp"
#include "wrapper/selector.hpp"

#include "cxxopts.hpp"

#ifdef QUALITY
#include "operation_log.hpp"
#else
#ifdef WITH_PAPI
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
#endif

using pq_type = PQWrapper<>;
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
    std::filesystem::path thread_stat_file;
#ifdef QUALITY
    std::filesystem::path log_file;
    std::filesystem::path metrics_file;
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
                    return "push";
                case Mode::PushAscending:
                    return "push (ascending)";
                default:
                    return "";
            }
        };

        out << "Threads: " << settings.num_threads << '\n'
            << "Prefill per thread: " << settings.prefill_per_thread << '\n'
            << "Prefill range: [" << static_cast<unsigned long>(settings.min_prefill) << ", "
            << static_cast<unsigned long>(settings.max_prefill) << "]\n"
            << "Iterations: " << settings.num_iterations << '\n'
            << "Mode: " << mode_name() << '\n';
        if (settings.mode == Mode::Update) {
            out << "Update range: [" << settings.min_update << ", " << settings.max_update << "]\n";
        } else if (settings.mode == Mode::PushRandom || settings.mode == Mode::Random) {
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
#ifdef QUALITY
        if (!settings.log_file.empty()) {
            out << "Log operations to: " << settings.log_file << '\n';
        }
        if (!settings.metrics_file.empty()) {
            out << "Write per element metrics to: " << settings.metrics_file << '\n';
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

struct Stats {
    std::chrono::high_resolution_clock::time_point start_time{};
    std::chrono::high_resolution_clock::time_point end_time{};
    long long num_iterations = 0;
    long long num_failed_pops = 0;
#ifdef WITH_PAPI
    bool papi_success = true;
    std::vector<long long> papi_counter{};
#endif
#ifdef MQ_COUNT_STATS
    multiqueue::Counters mq_stats{};
#endif

    template <typename It>
    static Stats accumulate(It begin, It end) {
        if (begin == end) {
            return {};
        }
        Stats accum_stats = *begin;
        ++begin;
        for (; begin != end; ++begin) {
            accum_stats.start_time = std::min(accum_stats.start_time, begin->start_time);
            accum_stats.end_time = std::max(accum_stats.end_time, begin->end_time);
#ifdef WITH_PAPI
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
        }
        return accum_stats;
    }

    static void write_header(Settings const& settings, std::ostream& out) {
        (void)settings;
        out << "time,pushes,pops,failed_pops";
#ifdef MQ_COUNT_STATS
        out << ",locked-push-pq,empty-pop-pqs,locked-pop-pq,stale-pop-pq";
#endif
#ifdef WITH_PAPI
        for (auto const& name : settings.papi_events) {
            out << "," << name;
        }
#endif
    }

    static void write(Settings const& settings, Stats const& stats, std::ostream& out) {
        out << (stats.end_time - stats.start_time).count() << ','
            << stats.num_iterations * Settings::pushes_per_iteration(settings.mode) << ','
            << stats.num_iterations * Settings::pops_per_iteration(settings.mode) << ',' << stats.num_failed_pops;
#ifdef MQ_COUNT_STATS
        out << ',' << stats.mq_stats.locked_push_pq << ',' << stats.mq_stats.empty_pop_pqs << ','
            << stats.mq_stats.locked_pop_pq << ',' << stats.mq_stats.stale_pop_pq;
#endif
#ifdef WITH_PAPI
        for (auto c : stats.papi_counter) {
            out << ',' << c;
        }
#endif
    }
};

struct BenchmarkData {
    std::vector<std::vector<key_type>> prefill;
    std::vector<long long> work;
    std::vector<Stats> stats;
#ifdef QUALITY
    std::vector<operation_log::OperationLog> op_logs;
#endif
    std::chrono::steady_clock::time_point start_time;
};

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
void work_loop(Settings const& settings, task::Control const& tc, Stats& stats, Work&& work) {
    static std::atomic<long long> counter{0};
    tc.synchronize();
#ifdef WITH_PAPI
    int event_set = PAPI_NULL;
    bool papi_started = false;
    if (!settings.papi_events.empty()) {
        stats.papi_counter.resize(settings.papi_events.size());
        try {
            event_set = start_papi(settings.papi_events);
            papi_started = true;
        } catch (std::runtime_error const& e) {
            tc.write(std::cerr) << "Error: " << e.what() << '\n';
            stats.papi_success = false;
        }
    }
#endif
    stats.start_time = std::chrono::high_resolution_clock::now();
    auto from = counter.fetch_add(settings.batch_size, std::memory_order_relaxed);
    for (; from < settings.num_iterations - settings.batch_size;
         from = counter.fetch_add(settings.batch_size, std::memory_order_relaxed)) {
        for (auto i = 0; i < settings.batch_size; ++i) {
            work(from + i);
        }
        stats.num_iterations += settings.batch_size;
        if (settings.timeout > std::chrono::seconds::zero()) {
            auto now = std::chrono::high_resolution_clock::now();
            if (now - stats.start_time > settings.timeout) {
                from = settings.num_iterations;
                break;
            }
        }
    }
    for (; from < settings.num_iterations; ++from) {
        work(from);
        ++stats.num_iterations;
    }
    stats.end_time = std::chrono::high_resolution_clock::now();
#ifdef WITH_PAPI
    if (papi_started) {
        if (int ret = PAPI_stop(event_set, stats.papi_counter.data()); ret != PAPI_OK) {
            tc.write(std::cerr) << "Error: Failed to stop performance counters\n";
            stats.papi_success = false;
        }
    }
#endif
    tc.synchronize();
}

std::ostream& log(std::ostream& out, std::chrono::steady_clock::time_point start) {
    out << '[' << std::fixed << std::setprecision(2)
        << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() << "s] ";
    return out;
}

void benchmark_thread(Settings const& settings, task::Control tc, pq_type& pq, BenchmarkData& benchmark_data) {
    auto handle = pq.get_handle();
    Stats stats{};
#ifdef QUALITY
    operation_log::OperationLog op_log{};
    op_log.pushes.reserve(static_cast<std::size_t>(
        settings.prefill_per_thread + settings.num_iterations * Settings::pushes_per_iteration(settings.mode)));
    op_log.pops.reserve(
        static_cast<std::size_t>(settings.num_iterations * Settings::pops_per_iteration(settings.mode)));
#endif

    auto pq_push = [&](key_type key, value_type value) {
        handle.push({key, value});
#ifdef QUALITY
        auto tick = static_cast<long long>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        op_log.pushes.push_back({tick, key, static_cast<std::size_t>(value)});
#endif
    };

    auto pq_pop = [&]() {
        while (true) {
#ifdef QUALITY
            auto tick = static_cast<long long>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
#endif
            auto retval = handle.try_pop();
            if (retval) {
#ifdef QUALITY
                op_log.pops.push_back({tick, static_cast<std::size_t>(retval->second)});
#endif
                return retval->first;
            }
        }
    };

    tc.once([&]() { log(std::clog, benchmark_data.start_time) << "Generating work...\n"; });
    tc.synchronize();
    generate_work(settings, tc.id(), benchmark_data.prefill[static_cast<std::size_t>(tc.id())], benchmark_data.work);
    tc.synchronize();
    tc.once([&]() { log(std::clog, benchmark_data.start_time) << "Prefilling...\n"; });
    tc.synchronize();
    for (auto i = 0LL; i < settings.prefill_per_thread; ++i) {
        pq_push(benchmark_data.prefill[static_cast<std::size_t>(tc.id())][static_cast<std::size_t>(i)],
                static_cast<value_type>(tc.id() * settings.prefill_per_thread + i));
    }
    tc.synchronize();
    tc.once([&]() { log(std::clog, benchmark_data.start_time) << "Running benchmark...\n"; });
#ifdef MQ_COUNT_STATS
    handle.reset_counters();
#endif
    switch (settings.mode) {
        case Settings::Mode::Update: {
            work_loop(settings, tc, stats, [&](auto i) {
                auto key = pq_pop();
                auto new_key = static_cast<key_type>(static_cast<long long>(key) +
                                                     benchmark_data.work[static_cast<std::size_t>(i)]);
                pq_push(new_key, static_cast<value_type>(settings.num_threads * settings.prefill_per_thread + i));
            });
            break;
        }
        case Settings::Mode::Random: {
            work_loop(settings, tc, stats, [&](auto i) {
                pq_pop();
                auto new_key = static_cast<key_type>(benchmark_data.work[static_cast<std::size_t>(i)]);
                pq_push(new_key, static_cast<value_type>(settings.num_threads * settings.prefill_per_thread + i));
            });
            break;
        }
        case Settings::Mode::PushRandom:
        case Settings::Mode::PushAscending: {
            work_loop(settings, tc, stats, [&](auto i) {
                auto new_key = static_cast<key_type>(benchmark_data.work[static_cast<std::size_t>(i)]);
                pq_push(new_key, static_cast<value_type>(settings.num_threads * settings.prefill_per_thread + i));
            });
            break;
        }
        case Settings::Mode::Pop: {
            work_loop(settings, tc, stats, [&](auto) { pq_pop(); });
            break;
        }
    }
#ifdef MQ_COUNT_STATS
    stats.mq_stats = data.handle.get_counters();
#endif
    benchmark_data.stats[static_cast<std::size_t>(tc.id())] = stats;
#ifdef QUALITY
    benchmark_data.op_logs[static_cast<std::size_t>(tc.id())] = std::move(op_log);
#endif
}

bool run_benchmark(Settings const& settings, cxxopts::ParseResult const& parse_result) {
    BenchmarkData benchmark_data{};
    benchmark_data.start_time = std::chrono::steady_clock::now();

    std::ofstream thread_stat_out;
    if (!settings.thread_stat_file.empty()) {
        thread_stat_out = std::ofstream(settings.thread_stat_file);
        if (!thread_stat_out) {
            std::cerr << "Error: Could not open file " << settings.thread_stat_file << " for writing" << std::endl;
            return false;
        }
    }
#ifdef QUALITY
    std::ofstream metrics_out;
    if (!settings.metrics_file.empty()) {
        metrics_out = std::ofstream(settings.metrics_file);
        if (!metrics_out) {
            std::cerr << "Error: Could not open file " << settings.metrics_file << " for writing" << std::endl;
            return false;
        }
    }
    std::ofstream log_out;
    if (!settings.log_file.empty()) {
        log_out = std::ofstream(settings.log_file);
        if (!log_out) {
            std::cerr << "Error: Could not open file " << settings.log_file << " for writing" << std::endl;
            return false;
        }
    }
#endif

    benchmark_data.prefill.resize(static_cast<std::size_t>(settings.num_threads));
    benchmark_data.work.resize(
        static_cast<std::size_t>(settings.num_iterations * Settings::pushes_per_iteration(settings.mode)));
    benchmark_data.stats.resize(static_cast<std::size_t>(settings.num_threads));
#ifdef QUALITY
    benchmark_data.op_logs.resize(static_cast<std::size_t>(settings.num_threads));
#endif
    auto max_capacity = benchmark_data.prefill.size() +
        (settings.mode == Settings::Mode::PushAscending || settings.mode == Settings::Mode::PushRandom
             ? benchmark_data.work.size()
             : 0);
    auto pq = create(settings.num_threads, 2 * max_capacity, parse_result);
    std::clog << "Priority queue: ";
    describe(pq, std::clog) << '\n' << '\n';

    task::Runner runner(affinity::NUMA{cores_per_numa_node, num_numa_nodes}, settings.num_threads,
                        [&](auto tc) { benchmark_thread(settings, tc, pq, benchmark_data); });
    runner.wait();

    if (thread_stat_out.is_open()) {
        Stats::write_header(settings, thread_stat_out);
        thread_stat_out << '\n';
        for (auto const& stats : benchmark_data.stats) {
            Stats::write(settings, stats, thread_stat_out);
            thread_stat_out << '\n';
        }
        thread_stat_out.close();
    }
    auto total_stats = Stats::accumulate(benchmark_data.stats.begin(), benchmark_data.stats.end());
#ifdef QUALITY
    log(std::clog, benchmark_data.start_time) << "Merging operation logs...\n";
    operation_log::OperationLog op_log;
    op_log.pushes.reserve(std::accumulate(benchmark_data.op_logs.begin(), benchmark_data.op_logs.end(), std::size_t{0},
                                          [](std::size_t size, auto const& l) { return size + l.pushes.size(); }));
    op_log.pops.reserve(std::accumulate(benchmark_data.op_logs.begin(), benchmark_data.op_logs.end(), std::size_t{0},
                                        [](std::size_t size, auto const& l) { return size + l.pops.size(); }));
    for (auto const& l : benchmark_data.op_logs) {
        op_log.pushes.insert(op_log.pushes.end(), l.pushes.begin(), l.pushes.end());
        op_log.pops.insert(op_log.pops.end(), l.pops.begin(), l.pops.end());
    }
    std::sort(op_log.pushes.begin(), op_log.pushes.end());
    std::sort(op_log.pops.begin(), op_log.pops.end());
    if (log_out.is_open()) {
        log(std::clog, benchmark_data.start_time) << "Writing operation log...\n";
        operation_log::write(op_log, log_out);
        log_out.close();
    }
    log(std::clog, benchmark_data.start_time) << "Replaying operations...\n";
    auto metrics = operation_log::replay(op_log);
    if (metrics_out.is_open()) {
        metrics_out << "rank_error,delay\n";
        for (auto const& m : metrics) {
            metrics_out << m.rank_error << ',' << m.delay << '\n';
        }
        metrics_out.close();
    }
    assert(!metrics.empty());
    auto summed_metrics =
        std::reduce(metrics.begin(), metrics.end(), operation_log::Metrics{}, [](auto a, auto const& b) {
            a.rank_error += b.rank_error;
            a.delay += b.delay;
            return a;
        });
#endif
    log(std::clog, benchmark_data.start_time) << "Finished\n";
    Stats::write_header(settings, std::cout);
#ifdef QUALITY
    std::cout << ",rank-error,delay";
#endif
    std::cout << '\n';
    Stats::write(settings, total_stats, std::cout);
#ifdef QUALITY
    std::cout << ',' << summed_metrics.rank_error << ',' << summed_metrics.delay;
#endif
    std::cout << '\n';
    return true;
}

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
        ("thread-stats", "File to write thread stats to", cxxopts::value<std::filesystem::path>(settings.thread_stat_file), "PATH")
#ifdef QUALITY
        ("l,log-file", "File to write the operation log to", cxxopts::value<std::filesystem::path>(settings.log_file), "PATH")
        ("metrics-file", "File to write single metrics to", cxxopts::value<std::filesystem::path>(settings.metrics_file), "PATH")
#endif
#ifdef WITH_PAPI
        ("r,pc", "Performance counters", cxxopts::value<std::vector<std::string>>(settings.papi_events))
#endif
        // clang-format on
        ;
    add_options(cmd);
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
        settings.mode = parse_mode(args["mode"].as<char>());
    }
    settings.timeout = std::chrono::seconds(timeout_ms);

    std::clog << '\n';
    std::clog << "Command line:";
    for (int i = 0; i < argc; ++i) {
        std::clog << ' ' << argv[i];
    }
    std::clog << '\n';
    Settings::write(settings, std::clog);

#ifdef WITH_PAPI
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

    bool success = run_benchmark(settings, args);
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
