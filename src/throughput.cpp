#include "priority_queue_factory.hpp"
#include "thread_coordination.hpp"
#ifdef PQ_MQ
#include "multiqueue/config.hpp"
#endif

#include "cxxopts.hpp"

#include <atomic>
#include <cassert>
#include <chrono>
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

#ifdef USE_PAPI
extern "C" {
#include <papi.h>
#include <pthread.h>
}
#endif

using key_type = unsigned long;
using value_type = unsigned long;

using PriorityQueue = util::priority_queue_type<key_type, value_type, true>;

static constexpr std::size_t DefaultPrefillPerThread = 1 << 20;
static constexpr std::size_t DefaultOperationsPerThread = 1 << 24;
static constexpr int DefaultNumThreads = 4;
static constexpr double DefaultPopProbability = 0.5;
static constexpr key_type DefaultMaxKey = key_type{1} << 30;

struct Settings {
    std::size_t prefill_per_thread = DefaultPrefillPerThread;
    std::size_t operations_per_thread = DefaultOperationsPerThread;
    int num_threads = DefaultNumThreads;
    double pop_prob = DefaultPopProbability;
    int seed = 1;
    key_type min_key = 1;
    key_type max_key = DefaultMaxKey;
    bool no_work = false;
#ifdef USE_PAPI
    int use_pc = -1;
    int l3_groupsize = 4;
#endif
};

Settings settings;

#ifdef USE_PAPI
std::array const event_names{"PAPI_L1_DCA",
                             "PAPI_L2_DCM",
                             "PAPI_L2_DCH",
                             "PAPI_L2_DCR",
                             "perf::PERF_COUNT_HW_CACHE_L1D",
                             "perf::L1-DCACHE-LOADS",
                             "perf::L1-DCACHE-LOAD-MISSES",
                             "perf::L1-DCACHE-STORES",
                             "perf::L1-DCACHE-STORE-MISSES",
                             "perf_raw:r0729",
                             "perf_raw:rc860",
                             "perf_raw:r0864",
                             "perf::LLC-LOADS",
                             "perf::LLC-LOAD-MISSES",
                             "perf::LLC-STORES",
                             "perf::LLC-STORE-MISSES"};
#endif

class Benchmark {
   public:
    using Handle = typename PriorityQueue::Handle;

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

    struct Data {
        std::vector<Operation> operations;
        thread_coordination::duration_type prefill_time{};
        thread_coordination::duration_type work_time{};
        std::atomic_llong num_failed_pops{0};
#ifdef USE_PAPI
        std::atomic_llong counter{0};
#endif
#ifdef MQ_COUNT_STATS
        std::atomic_size_t num_locking_failed{0};
        std::atomic_size_t num_resets{0};
        std::atomic_size_t use_counts{0};
#endif
        explicit Data(std::size_t num_operations) : operations(num_operations) {
        }
    };

   private:
#ifdef USE_PAPI
    static bool start_performance_counter(int& event_set) {
        if (int ret = PAPI_register_thread(); ret != PAPI_OK) {
            return false;
        }
        if (int ret = PAPI_create_eventset(&event_set); ret != PAPI_OK) {
            return false;
        }
        int code{};
        if (int ret = PAPI_event_name_to_code(event_names[static_cast<std::size_t>(settings.use_pc)], &code);
            ret != PAPI_OK) {
            return false;
        }
        if (int ret = PAPI_add_event(event_set, code); ret != PAPI_OK) {
            return false;
        }
        if (int ret = PAPI_start(event_set); ret != PAPI_OK) {
            return false;
        }
        return true;
    }
#endif

    static void work(thread_coordination::Context const& ctx, Handle& handle, Data& data) {
        long long num_failed_pops = 0;
#ifdef USE_PAPI
        bool papi_started = false;
        int event_set = PAPI_NULL;
        if (settings.use_pc >= 0 &&
            (static_cast<std::size_t>(settings.use_pc) < event_names.size() - 4 ||
             ctx.get_id() % settings.l3_groupsize == 0)) {
            papi_started = start_performance_counter(event_set);
            if (!papi_started) {
                ctx.write(std::cerr) << "Failed to start counters\n";
            }
        }
#endif
        volatile PriorityQueue::mapped_type val;
        PriorityQueue::value_type retval;
        ctx.execute_synchronized_blockwise_timed(
            data.work_time, data.operations.begin(), data.operations.end(),
            [&handle, &num_failed_pops, &val, &retval](auto block_begin, auto block_end) {
                for (auto it = block_begin; it != block_end; ++it) {
                    if (it->is_pop()) {
                        while (!handle.try_pop(retval)) {
                            ++num_failed_pops;
                        }
                        val = retval.second;
                    } else {
                        key_type key = it->get_insert_key();
                        handle.push({key, key});
                    }
                }
            });
#ifdef USE_PAPI
        if (papi_started) {
            long long counter{};
            if (int ret = PAPI_stop(event_set, &counter); ret != PAPI_OK) {
                ctx.write(std::cerr) << "Failed to stop counters\n";
            } else {
                data.counter += counter;
            }
        }
#endif
        data.num_failed_pops += num_failed_pops;
    }

   public:
    static void run(thread_coordination::Context ctx, Data& data, PriorityQueue& pq) {
        if (ctx.get_id() == 0) {
            std::clog << "Generating operations..." << std::flush;
        }

        std::seed_seq seed{settings.seed, ctx.get_id()};
        std::default_random_engine rng(seed);
        auto key_dist = std::uniform_int_distribution<key_type>(settings.min_key, settings.max_key);
        auto pop_dist = std::bernoulli_distribution(settings.pop_prob);

        // Prefill
        std::vector<key_type> prefill_keys(settings.prefill_per_thread);
        std::generate(prefill_keys.begin(), prefill_keys.end(), [&key_dist, &rng]() { return key_dist(rng); });

        // Workload
        auto block_begin =
            data.operations.begin() + ctx.get_id() * static_cast<std::ptrdiff_t>(settings.operations_per_thread);
        std::generate_n(block_begin, settings.operations_per_thread,
                        [&]() { return pop_dist(rng) ? Operation{} : Operation{key_dist(rng)}; });

        ctx.synchronize([]() { std::clog << "done\nPrefilling..." << std::flush; });

        Handle handle = pq.get_handle(ctx.get_id());

        ctx.execute_synchronized_timed(data.prefill_time, [&handle, &prefill_keys]() {
            for (auto k : prefill_keys) {
                handle.push({k, k});
            }
        });

        if (ctx.get_id() == 0) {
            std::clog << "done\nWorking..." << std::flush;
        }
        if (!settings.no_work) {
#ifdef MQ_COUNT_STATS
            handle.stats.num_locking_failed = 0;
            handle.stats.num_resets = 0;
            handle.stats.use_counts = 0;
#endif
            work(ctx, handle, data);
#ifdef MQ_COUNT_STATS
            data.num_locking_failed += handle.stats.num_locking_failed;
            data.num_resets += handle.stats.num_resets;
            data.use_counts += handle.stats.use_counts;
#endif
        }
        if (ctx.get_id() == 0) {
            std::clog << "done\n" << std::endl;
        }
    }
};

int main(int argc, char* argv[]) {
#ifndef NDEBUG
    std::clog << "Build type: Debug\n";
#else
    std::clog << "Build type: Release\n";
#endif
#ifdef USE_PAPI
    std::clog << "Performance counter: enabled\n";
#else
    std::clog << "Performance counter: disabled\n";
#endif
    std::clog << "L1 cache linesize (bytes): " << L1_CACHE_LINESIZE << '\n';
    std::clog << "Pagesize (bytes): " << PAGESIZE << '\n';
    std::clog << '\n';

#ifdef PQ_MQ
    multiqueue::Config mq_config;
#endif

    cxxopts::Options options(argv[0]);
    options.add_options()
        // clang-format off
        ("h,help", "Print this help")
        ("j,threads", "The number of threads", cxxopts::value<int>(settings.num_threads), "NUMBER")
        ("p,prefill", "The prefill per thread", cxxopts::value<std::size_t>(settings.prefill_per_thread), "NUMBER")
        ("n,ops", "The number of operations per thread", cxxopts::value<std::size_t>(settings.operations_per_thread), "NUMBER")
        ("d,pop-prob", "Specify the probability of pops in [0,1]", cxxopts::value<double>(settings.pop_prob), "NUMBER")
        ("l,min", "Specify the min key", cxxopts::value<key_type>(settings.min_key), "NUMBER")
        ("m,max", "Specify the max key", cxxopts::value<key_type>(settings.max_key), "NUMBER")
        ("s,seed", "Specify the initial seed", cxxopts::value<int>(settings.seed), "NUMBER")
#ifdef PQ_MQ
        ("c,factor", "The factor for queues", cxxopts::value<int>(mq_config.c), "NUMBER")
        ("k,stickiness", "The stickiness period", cxxopts::value<int>(mq_config.stickiness), "NUMBER")
#endif
        ("x,no-work", "Don't perform the actual benchmark", cxxopts::value<bool>(settings.no_work))
#ifdef USE_PAPI
        ("r,pc", "Performance counter to count", cxxopts::value<int>(settings.use_pc))
#endif
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

    std::clog << "Command line:";
    for (int i = 0; i < argc; ++i) {
        std::clog << ' ' << argv[i];
    }
    std::clog << '\n';
    std::clog << "Threads: " << settings.num_threads << '\n'
              << "Prefill per thread: " << settings.prefill_per_thread << '\n'
              << "Operations per thread: " << settings.operations_per_thread << '\n'
              << "Pop probability: " << std::fixed << std::setprecision(2) << settings.pop_prob << '\n'
              << "Min key: " << settings.min_key << '\n'
              << "Max key: " << settings.max_key << '\n'
              << "Seed: " << settings.seed << '\n'
#ifdef PQ_MQ
              << "Factor: " << mq_config.c << '\n'
              << "Stickiness: " << mq_config.stickiness << '\n'
#endif
              << '\n';

#ifdef USE_PAPI
    if (settings.use_pc >= 0) {
        if (static_cast<std::size_t>(settings.use_pc) >= event_names.size()) {
            std::cerr << "Invalid event!\n";
            return 1;
        }
        if (int ret = PAPI_library_init(PAPI_VER_CURRENT); ret != PAPI_VER_CURRENT) {
            std::cerr << "Error initializing PAPI\n";
            return 1;
        }
        if (int ret = PAPI_thread_init((unsigned long (*)())(pthread_self)); ret != PAPI_OK) {
            std::cerr << "Error initializing threads for PAPI\n";
            return 1;
        }
        if (int ret = PAPI_query_named_event(event_names[static_cast<std::size_t>(settings.use_pc)]); ret != PAPI_OK) {
            std::cerr << "Cannot measure event '" << event_names[static_cast<std::size_t>(settings.use_pc)] << "'\n";
            return 1;
        }
    }
#endif

#ifdef PQ_MQ
    mq_config.seed = settings.seed;
    auto pq = PriorityQueue(settings.num_threads,
                            settings.prefill_per_thread * static_cast<unsigned int>(settings.num_threads), mq_config);
#else
    auto pq = PriorityQueue(settings.num_threads);
#endif

    std::clog << "Data structure: ";
    util::describe::describe(std::clog, pq);
    std::clog << '\n';

    Benchmark::Data benchmark_data(static_cast<std::size_t>(settings.num_threads) * settings.operations_per_thread);

    thread_coordination::TaskHandle task_handle(settings.num_threads, Benchmark::run, std::ref(benchmark_data),
                                                std::ref(pq));

    task_handle.wait();

    std::clog << "Prefill time (s): " << std::setprecision(3)
              << std::chrono::duration<double>(benchmark_data.prefill_time).count() << '\n';
    std::clog << "Work time (s): " << std::setprecision(3)
              << std::chrono::duration<double>(benchmark_data.work_time).count() << '\n';
    std::clog << "Failed pops: " << benchmark_data.num_failed_pops << '\n';
#ifdef USE_PAPI
    if (settings.use_pc >= 0) {
        std::clog << event_names[static_cast<std::size_t>(settings.use_pc)] << ": " << benchmark_data.counter << '\n';
    }
#endif
#ifdef MQ_COUNT_STATS
    std::clog << "Failed locks per operation: "
              << static_cast<double>(benchmark_data.num_locking_failed) /
            static_cast<double>(settings.num_threads * settings.operations_per_thread)
              << '\n';
    std::clog << "Average queue use count: "
              << static_cast<double>(benchmark_data.use_counts) / static_cast<double>(benchmark_data.num_resets)
              << '\n';
#endif

    std::cout << settings.num_threads << ',' << settings.prefill_per_thread << ',' << settings.operations_per_thread
              << ',' << settings.pop_prob << ',' << settings.min_key << ',' << settings.max_key << ',' << settings.seed
              << ',' << std::fixed << std::setprecision(3)
              << std::chrono::duration<double>(benchmark_data.prefill_time).count() << ',' << std::fixed
              << std::setprecision(3) << std::chrono::duration<double>(benchmark_data.work_time).count() << ','
              << benchmark_data.num_failed_pops;
#ifdef USE_PAPI
    if (settings.use_pc >= 0) {
        std::cout << ',' << benchmark_data.counter;
    } else {
        std::cout << ",n/a";
    }
#else
    std::cout << ",n/a";
#endif
#ifdef MQ_COUNT_STATS
    std::cout << ',' << benchmark_data.num_resets << ',' << benchmark_data.use_counts;
#else
    std::cout << ",n/a,n/a";
#endif
    std::cout << std::endl;
    return 0;
}
