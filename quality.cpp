#include "common/priority_queue_factory.hpp"
#include "common/thread_coordination.hpp"
#ifdef PQ_MQ
#include "multiqueue/config.hpp"
#endif

#include "cxxopts.hpp"

#include <x86intrin.h>
#include <atomic>
#include <cassert>
#include <chrono>
#include <ctime>
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
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

using key_type = unsigned long;
using value_type = unsigned long;

using PriorityQueue = util::priority_queue_type<key_type, value_type, true>;

static constexpr unsigned int bits_for_value = std::numeric_limits<value_type>::digits - 8;
static constexpr unsigned int MaxNumThreads = 1UL << (std::numeric_limits<value_type>::digits - bits_for_value);
static constexpr value_type value_mask = (value_type{1} << bits_for_value) - 1;

static value_type to_value(unsigned int thread_id, value_type elem_id) noexcept {
    assert(static_cast<value_type>(thread_id) <
           (value_type{1} << (std::numeric_limits<value_type>::digits - bits_for_value)));
    assert(elem_id < (value_type{1} << bits_for_value));
    return (value_type{thread_id} << (bits_for_value)) | (elem_id & value_mask);
}

static unsigned int get_thread_id(value_type value) noexcept {
    return static_cast<unsigned int>(value >> bits_for_value);
}

static value_type get_elem_id(value_type value) noexcept {
    return value & value_mask;
}

static constexpr std::size_t DefaultPrefill = 1 << 20;
static constexpr std::size_t DefaultOperations = 1 << 24;
static constexpr unsigned int DefaultNumThreads = 4;
static constexpr double DefaultPopProbability = 0.5;
static constexpr key_type DefaultMaxKey = key_type{1} << 30;

struct Settings {
    std::size_t prefill = DefaultPrefill;
    std::size_t operations = DefaultOperations;
    unsigned int num_threads = DefaultNumThreads;
    double pop_prob = DefaultPopProbability;
    std::uint32_t seed = 1;
    key_type min_key = 1;
    key_type max_key = DefaultMaxKey;
    std::chrono::microseconds sleep_between_operations = std::chrono::microseconds::zero();
};

class Benchmark {
   public:
    using Handle = typename PriorityQueue::handle_type;

    using tick_type = std::uint64_t;

    class Operation {
        key_type data;

       public:
        static constexpr key_type PopOp = 0;

        Operation() : data{PopOp} {
        }

        Operation(key_type insert_key) : data{insert_key} {
            assert(insert_key != PopOp);
        }

        [[nodiscard]] constexpr bool is_pop() const noexcept {
            return data == PopOp;
        }
        [[nodiscard]] constexpr key_type get_insert_key() const noexcept {
            assert(!is_pop());
            return data;
        }
    };

    struct PushLogEntry {
        tick_type tick;
        key_type key;
    };

    struct PopLogEntry {
        tick_type tick;
        value_type value;
    };

    struct Data {
        std::vector<Operation> operations;
        thread_coordination::duration_type prefill_time{};
        thread_coordination::duration_type work_time{};
        std::vector<std::vector<PushLogEntry>> push_log;
        std::vector<std::vector<PopLogEntry>> pop_log;
        std::vector<std::vector<tick_type>> failed_pop_log;
        std::atomic_size_t cache_accesses = 0;
        std::atomic_size_t cache_misses = 0;
        std::mutex m;
    };

   private:
    static tick_type get_tick() noexcept {
#ifdef USE_TSC
        // Not synchronized among sockets
        _mm_lfence();
        return __rdtsc();
        _mm_lfence();
#else
        static constexpr long nsec_per_s = 1000000000;
        timespec ts;
        clock_gettime(CLOCK_REALTIME, &ts);
        return static_cast<tick_type>(ts.tv_sec * nsec_per_s + ts.tv_nsec);
#endif
    }

    static void work(thread_coordination::Context const& ctx, Handle& handle, Data& data) {
        ctx.execute_synchronized_blockwise_timed(
            data.work_time, data.operations.begin(), data.operations.end(), [&](auto block_begin, auto block_end) {
                for (auto it = block_begin; it != block_end; ++it) {
                    if (it->is_pop()) {
                        PriorityQueue::value_type retval;
                        bool success = handle.try_pop(retval);
                        auto tick = get_tick();
                        if (success) {
                            data.pop_log[ctx.get_id()].push_back({tick, retval.second});
                        } else {
                            data.failed_pop_log[ctx.get_id()].push_back(tick);
                        }
                    } else {
                        key_type key = it->get_insert_key();
                        value_type value = to_value(ctx.get_id(), data.push_log[ctx.get_id()].size());
                        handle.push({key, value});
                        auto tick = get_tick();
                        data.push_log[ctx.get_id()].push_back({tick, key});
                    }
                }
            });
    }

    static void work(thread_coordination::Context const& ctx, Handle& handle, Data& data,
                     std::chrono::microseconds sleep, std::default_random_engine& rng) {
        auto sleep_dist = std::uniform_int_distribution<long>(1, sleep.count());
        ctx.execute_synchronized_blockwise_timed(
            data.work_time, data.operations.begin(), data.operations.end(), [&](auto block_begin, auto block_end) {
                for (auto it = block_begin; it != block_end; ++it) {
                    if (it->is_pop()) {
                        PriorityQueue::value_type retval;
                        bool success = handle.try_pop(retval);
                        auto tick = get_tick();
                        if (success) {
                            data.pop_log[ctx.get_id()].push_back({tick, retval.second});
                        } else {
                            data.failed_pop_log[ctx.get_id()].push_back(tick);
                        }
                    } else {
                        key_type key = it->get_insert_key();
                        value_type value = to_value(ctx.get_id(), data.push_log[ctx.get_id()].size());
                        handle.push({key, value});
                        auto tick = get_tick();
                        data.push_log[ctx.get_id()].push_back({tick, key});
                    }
                    std::this_thread::sleep_for(std::chrono::microseconds{sleep_dist(rng)});
                }
            });
    }

   public:
    static void run(thread_coordination::Context ctx, Settings const& settings, Data& data, PriorityQueue& pq) {
        if (ctx.get_id() == 0) {
            std::clog << "Generating operations..." << std::flush;
        }

        std::seed_seq seed{settings.seed, static_cast<std::uint32_t>(ctx.get_id())};
        std::default_random_engine rng(seed);
        auto key_dist = std::uniform_int_distribution<key_type>(settings.min_key, settings.max_key);
        auto pop_dist = std::bernoulli_distribution(settings.pop_prob);

        // Prefill
        std::vector<key_type> prefill_keys(settings.prefill / settings.num_threads);
        std::generate(prefill_keys.begin(), prefill_keys.end(), [&key_dist, &rng]() { return key_dist(rng); });

        // Workload
        ctx.execute_synchronized_blockwise(data.operations.begin(), data.operations.end(),
                                           [&](auto block_begin, auto block_end) {
                                               for (auto it = block_begin; it != block_end; ++it) {
                                                   *it = pop_dist(rng) ? Operation{} : Operation{key_dist(rng)};
                                               }
                                           });
        if (ctx.get_id() == 0) {
            std::clog << "done\nPrefilling..." << std::flush;
        }

        data.pop_log[ctx.get_id()].reserve(settings.operations);
        data.push_log[ctx.get_id()].reserve(settings.prefill);
        Handle handle = [&data, &pq]() {
            std::scoped_lock l{data.m};
            return pq.get_handle();
        }();

        ctx.execute_synchronized_timed(data.prefill_time, [id = ctx.get_id(), &handle, &prefill_keys, &data]() {
            for (auto k : prefill_keys) {
                value_type value = to_value(id, data.push_log[id].size());
                handle.push({k, value});
                data.push_log[id].push_back({0, k});
            }
        });

        if (ctx.get_id() == 0) {
            std::clog << "done\nWorking..." << std::flush;
        }

        if (settings.sleep_between_operations == std::chrono::microseconds::zero()) {
            work(ctx, handle, data);
        } else {
            work(ctx, handle, data, settings.sleep_between_operations, rng);
        }

        if (ctx.get_id() == 0) {
            std::clog << "done\n" << std::endl;
        }
    }
};

int main(int argc, char* argv[]) {
#ifndef NDEBUG
    std::cout << "Build type: Debug\n";
#else
    std::cout << "Build type: Release\n";
#endif
#ifdef USE_PAPI
    std::cout << "Performance counter: enabled\n";
#else
    std::cout << "Performance counter: disabled\n";
#endif
    std::cout << "L1 cache linesize (bytes): " << L1_CACHE_LINESIZE << '\n';
    std::cout << "Pagesize (bytes): " << PAGESIZE << '\n';
    std::cout << '\n';

    Settings settings;
#ifdef PQ_MQ
    multiqueue::Config mq_config;
#endif

    std::filesystem::path out_file;
    std::filesystem::path log_file;

    cxxopts::Options options(argv[0]);
    options.add_options()
        // clang-format off
        ("h,help", "Print this help")
        ("j,threads", "The number of threads", cxxopts::value<unsigned int>(settings.num_threads), "NUMBER")
        ("p,prefill", "The prefill", cxxopts::value<std::size_t>(settings.prefill), "NUMBER")
        ("n,ops", "The number of operations", cxxopts::value<std::size_t>(settings.operations), "NUMBER")
        ("d,pop-prob", "Specify the probability of pops in [0,1]", cxxopts::value<double>(settings.pop_prob), "NUMBER")
        ("l,min", "Specify the min key", cxxopts::value<key_type>(settings.min_key), "NUMBER")
        ("m,max", "Specify the max key", cxxopts::value<key_type>(settings.max_key), "NUMBER")
        ("s,seed", "Specify the initial seed", cxxopts::value<std::uint32_t>(settings.seed), "NUMBER")
#ifdef PQ_MQ
        ("c,factor", "The factor for queues", cxxopts::value<unsigned int>(mq_config.c), "NUMBER")
        ("k,stickiness", "The stickiness period", cxxopts::value<unsigned int>(mq_config.stickiness), "NUMBER")
#endif
        ("o,outfile", "Output data in csv (comma-separated)", cxxopts::value<std::filesystem::path>(out_file), "PATH")
        ("f,log-file", "Path to write the log file to", cxxopts::value<std::filesystem::path>(log_file)->default_value("log.txt"), "PATH")
        // clang-format on
        ;

    try {
        auto result = options.parse(argc, argv);
        if (result.count("help") > 0) {
            std::cerr << options.help() << std::endl;
            return 0;
        }
    } catch (cxxopts::OptionParseException const& e) {
        std::cerr << "Error parsing arguments: " << e.what() << '\n';
        std::cerr << options.help() << std::endl;
        return 1;
    }

    std::cout << "Command line:";
    for (int i = 0; i < argc; ++i) {
        std::cout << ' ' << argv[i];
    }
    std::cout << '\n';
    std::cout << "Threads: " << settings.num_threads << '\n'
              << "Prefill: " << settings.prefill << '\n'
              << "Operations: " << settings.operations << '\n'
              << "Pop probability: " << std::fixed << std::setprecision(2) << settings.pop_prob << '\n'
              << "Min key: " << settings.min_key << '\n'
              << "Max key: " << settings.max_key << '\n'
              << "Seed: " << settings.seed << '\n'
#ifdef PQ_MQ
              << "Factor: " << mq_config.c << '\n'
              << "Stickiness: " << mq_config.stickiness << '\n'
#endif
              << "Sleep between operations: " << settings.sleep_between_operations.count() << " us\n"
              << '\n';

    if (settings.num_threads > MaxNumThreads) {
        std::cerr << "Too many threads!" << std::endl;
        return 1;
    }

#ifdef PQ_MQ
    mq_config.seed = settings.seed;
    auto pq = PriorityQueue(settings.num_threads, mq_config);
#else
    auto pq = PriorityQueue(settings.num_threads);
#endif

    std::cout << "Data structure: ";
    util::describe::describe(std::cout, pq);
    std::cout << '\n';

    Benchmark::Data benchmark_data{};
    benchmark_data.operations.resize(settings.operations);
    benchmark_data.pop_log.resize(settings.num_threads);
    benchmark_data.push_log.resize(settings.num_threads);
    benchmark_data.failed_pop_log.resize(settings.num_threads);

    thread_coordination::TaskHandle<Benchmark> task_handle(settings.num_threads, std::cref(settings),
                                                           std::ref(benchmark_data), std::ref(pq));

    task_handle.join();

    std::cout << "Prefill time (s): " << std::setprecision(3)
              << std::chrono::duration<double>(benchmark_data.prefill_time).count() << '\n';
    std::cout << "Work time (s): " << std::setprecision(3)
              << std::chrono::duration<double>(benchmark_data.work_time).count() << '\n';

    std::clog << "Writing log file..." << std::flush;
    auto log_out = std::ofstream{log_file};
    if (!log_out.is_open()) {
        std::cerr << "Could not open log file\n";
    } else {
        log_out << "# " << settings.num_threads << ' ' << settings.prefill << '\n';
        for (unsigned int t = 0; t < settings.num_threads; ++t) {
            for (auto const& [tick, key] : benchmark_data.push_log[t]) {
                log_out << "i " << t << ' ' << tick << ' ' << key << '\n';
            }
        }

        for (unsigned int t = 0; t < settings.num_threads; ++t) {
            for (auto const& [tick, value] : benchmark_data.pop_log[t]) {
                log_out << "d " << t << ' ' << tick << ' ' << get_thread_id(value) << ' ' << get_elem_id(value) << '\n';
            }
        }

        for (unsigned int t = 0; t < settings.num_threads; ++t) {
            for (auto const& tick : benchmark_data.failed_pop_log[t]) {
                log_out << "f " << t << ' ' << tick << '\n';
            }
        }

        log_out << std::flush;
        log_out.close();
        std::clog << "done" << std::endl;
    }
    if (!out_file.empty()) {
        auto out = std::ofstream{out_file};
        if (!out.is_open()) {
            std::cerr << "Could not open file to write out benchmark_datas\n";
            return 1;
        }
        out << "threads,prefill,operations,pop_prob,min_key,max_key,"
               "seed,prefill_time,work_time\n";
        out << settings.num_threads << ',' << settings.prefill << ',' << settings.operations << ',' << settings.pop_prob
            << ',' << settings.min_key << ',' << settings.max_key << ',' << settings.seed << ',' << std::fixed
            << std::setprecision(3) << std::chrono::duration<double>(benchmark_data.prefill_time).count() << ','
            << std::fixed << std::setprecision(3) << std::chrono::duration<double>(benchmark_data.work_time).count()
            << ',' << std::endl;
    }
    return 0;
}
