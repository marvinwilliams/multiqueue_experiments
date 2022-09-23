#include <mutex>
#include "utils/priority_queue_factory.hpp"
#include "utils/thread_coordination.hpp"
#ifdef PQ_MQ
#include "multiqueue/config.hpp"
#endif

#include "external/cxxopts.hpp"

#ifdef USE_PAPI
extern "C" {
#include <papi.h>
#include <pthread.h>
}
#endif

#include <atomic>
#include <cassert>
#include <chrono>
#include <condition_variable>
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

using key_type = unsigned long;
using value_type = unsigned long;

class Operation {
    key_type data;

   public:
    static constexpr key_type PopOp = 0;

    explicit Operation() : data{PopOp} {
    }

    explicit Operation(key_type insert_key) : data{insert_key} {
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

using PriorityQueue = util::priority_queue_type<key_type, value_type, true>;
using Handle = typename PriorityQueue::handle_type;

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
    std::uint32_t seed = 1;
    key_type min_key = 1;
    key_type max_key = DefaultMaxKey;
};

class Benchmark {
   public:
    struct result_type {
        bool fail = false;
        thread_coordination::duration_type prefill_time{};
        thread_coordination::duration_type work_time{};
        std::size_t num_failed_pops = 0;
#ifdef USE_PAPI
        std::size_t cache_accesses = 0;
        std::size_t cache_misses = 0;
#endif
    };

    struct Data {
        bool fail = false;
        std::vector<Operation> operations;
        thread_coordination::duration_type prefill_time{};
        thread_coordination::duration_type work_time{};
        std::atomic_size_t num_failed_pops = 0;
#ifdef USE_PAPI
        std::atomic_size_t cache_accesses = 0;
        std::atomic_size_t cache_misses = 0;
#endif
        std::mutex m;
    };

   private:
    static bool check_operations(std::size_t prefill_elements, std::vector<Operation>& ops) {
        for (auto op : ops) {
            if (op.is_pop()) {
                if (prefill_elements == 0) {
                    return false;
                }
                --prefill_elements;
            } else {
                ++prefill_elements;
            }
        }
        return true;
    }

#ifdef USE_PAPI
    static bool start_performance_counter(int& event_set) {
        if (int ret = PAPI_register_thread(); ret != PAPI_OK) {
            return false;
        }
        if (int ret = PAPI_create_eventset(&event_set); ret != PAPI_OK) {
            return false;
        }
        if (int ret = PAPI_add_event(event_set, PAPI_L1_DCA); ret != PAPI_OK) {
            return false;
        }
        if (int ret = PAPI_add_named_event(event_set, "perf::L1-DCACHE-LOAD-MISSES"); ret != PAPI_OK) {
            return false;
        }
        if (int ret = PAPI_start(event_set); ret != PAPI_OK) {
            return false;
        }
        return true;
    }
#endif

    static void work(thread_coordination::Context const& ctx, Handle& handle, Data& data) {
        std::size_t num_failed_pops = 0;
#ifdef USE_PAPI
        int event_set = PAPI_NULL;
        bool papi_started = start_performance_counter(event_set);
        if (!papi_started) {
            ctx.write(std::cerr) << "Failed to start counters\n";
        }
#endif
        ctx.execute_synchronized_blockwise_timed(data.work_time, data.operations.begin(), data.operations.end(),
                                                 [&handle, &num_failed_pops](auto block_begin, auto block_end) {
                                                     for (auto it = block_begin; it != block_end; ++it) {
                                                         if (it->is_pop()) {
                                                             PriorityQueue::value_type retval;
                                                             // Let retval escape
                                                             asm volatile("" ::"g"(&retval));
                                                             while (!handle.try_pop(retval)) {
                                                                 ++num_failed_pops;
                                                             }
                                                             // "Use" memory to force write to retval
                                                             asm volatile("" ::: "memory");
                                                         } else {
                                                             key_type key = it->get_insert_key();
                                                             handle.push({key, key});
                                                         }
                                                     }
                                                 });
#ifdef USE_PAPI
        if (papi_started) {
            long long cache_stats[2];
            if (int ret = PAPI_stop(event_set, cache_stats); ret != PAPI_OK) {
                ctx.write(std::cerr) << "Failed to stop counters\n";
            }
            data.cache_accesses += cache_stats[0];
            data.cache_misses += cache_stats[1];
        }
#endif
        data.num_failed_pops += num_failed_pops;
    }

   public:
    static void run(thread_coordination::Context ctx, std::promise<result_type>& promise, Settings const& settings,
                    Data& data, PriorityQueue& pq) {
        std::seed_seq seed{settings.seed, static_cast<std::uint32_t>(ctx.get_id())};
        std::default_random_engine rng(seed);
        auto key_dist = std::uniform_int_distribution<key_type>(settings.min_key, settings.max_key);
        auto pop_dist = std::bernoulli_distribution(settings.pop_prob);

        if (ctx.get_id() == 0) {
            std::clog << "Generating operations..." << std::flush;
        }

        // Prefill
        std::vector<key_type> prefill_keys(settings.prefill_per_thread);
        std::generate(prefill_keys.begin(), prefill_keys.end(), [&key_dist, &rng]() { return key_dist(rng); });

        // Workload
        auto block_begin =
            data.operations.begin() + ctx.get_id() * static_cast<std::ptrdiff_t>(settings.operations_per_thread);
        std::generate_n(block_begin, settings.operations_per_thread,
                        [&]() { return pop_dist(rng) ? Operation{} : Operation{key_dist(rng)}; });

        ctx.synchronize([&settings, &data, &promise]() {
            std::clog << "done\nChecking operations..." << std::flush;
            std::size_t num_prefill_elements =
                static_cast<std::size_t>(settings.num_threads) * settings.prefill_per_thread;
            if (!check_operations(num_prefill_elements, data.operations)) {
                std::cerr << "\nError: Invalid operation sequence, increase prefill or decrease pop probability\n";
                data.fail = true;
                result_type result;
                result.fail = true;
                promise.set_value(result);

            } else {
                std::clog << "done\nPrefilling..." << std::flush;
            }
        });

        ctx.synchronize();

        if (data.fail) {
            return;
        }

        Handle handle = [&data, &pq]() {
            std::scoped_lock l{data.m};
            return pq.get_handle();
        }();

        ctx.execute_synchronized_timed(data.prefill_time, [&handle, &prefill_keys]() {
            for (auto k : prefill_keys) {
                handle.push({k, k});
            }
        });

        if (ctx.get_id() == 0) {
            std::clog << "done\nWorking..." << std::flush;
        }

        work(ctx, handle, data);

        if (ctx.get_id() == 0) {
            std::clog << "done\n" << std::endl;
            promise.set_value(result_type{data.fail, data.prefill_time, data.work_time, data.num_failed_pops
#ifdef USE_PAPI
                                          ,
                                          data.cache_accesses, data.cache_misses
#endif
            });
        }
    }
};

int main(int argc, char* argv[]) {
    std::cout << "Build config\n";
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

    std::cout << "Command line:";
    for (int i = 0; i < argc; ++i) {
        std::cout << ' ' << argv[i];
    }
    std::cout << "\n\n";

    Settings settings;
#ifdef PQ_MQ
    multiqueue::Config mq_config;
#endif
    std::filesystem::path out_file;

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
        ("s,seed", "Specify the initial seed", cxxopts::value<std::uint32_t>(settings.seed), "NUMBER")
#ifdef PQ_MQ
        ("c,factor", "The factor for queues", cxxopts::value<int>(mq_config.c), "NUMBER")
        ("k,stickiness", "The stickiness period", cxxopts::value<int>(mq_config.stickiness), "NUMBER")
#endif
        ("o,outfile", "Write the benchmark results to file", cxxopts::value<std::filesystem::path>(out_file), "PATH")
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

    std::cout << "Settings\n"
              << "Threads: " << settings.num_threads << '\n'
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
    if (int ret = PAPI_library_init(PAPI_VER_CURRENT); ret != PAPI_VER_CURRENT) {
        std::cerr << "Error initializing PAPI\n";
        return 1;
    }
    if (int ret = PAPI_thread_init((unsigned long (*)())(pthread_self)); ret != PAPI_OK) {
        std::cerr << "Error initializing threads for PAPI\n";
        return 1;
    }
    if (int ret = PAPI_query_event(PAPI_L1_DCA); ret != PAPI_OK) {
        std::cerr << "Cannot measure cache accesses\n";
        return 1;
    }
    if (int ret = PAPI_query_named_event("perf::L1-DCACHE-LOAD-MISSES"); ret != PAPI_OK) {
        std::cerr << "Cannot measure cache misses\n";
        return 1;
    }
#endif

#ifdef PQ_MQ
    mq_config.seed = settings.seed;
    auto pq = PriorityQueue(settings.num_threads, mq_config);
#else
    auto pq = PriorityQueue(settings.num_threads);
#endif

    std::cout << "Priority queue\nName: ";
    util::describe(std::cout, pq);
    std::cout << '\n';

    Benchmark::Data benchmark_data{};
    benchmark_data.operations.resize(static_cast<std::size_t>(settings.num_threads) * settings.operations_per_thread);

    thread_coordination::TaskHandle<Benchmark> task_handle(settings.num_threads);

    auto future = task_handle.run(thread_coordination::affinity::individual_cores{}, std::cref(settings),
                                  std::ref(benchmark_data), std::ref(pq));
    auto result = future.get();
    task_handle.join();
    if (result.fail) {
        std::cerr << "Benchmark failed\n";
        return 1;
    }

    std::cout << "Failed pops: " << result.num_failed_pops << '\n';
    std::cout << "Prefill time (s): " << std::setprecision(3)
              << std::chrono::duration<double>(result.prefill_time).count() << '\n';
    std::cout << "Work time (s): " << std::setprecision(3) << std::chrono::duration<double>(result.work_time).count()
              << '\n';
#ifdef USE_PAPI
    std::cout << "Cache accesses/misses: " << result.cache_accesses << '/' << result.cache_misses << '\n';
#else
    std::cout << "Cache accesses/misses: not measured\n";
#endif

    if (!out_file.empty()) {
        auto out = std::ofstream{out_file};
        if (!out.is_open()) {
            std::cerr << "Could not open file to write out results\n";
            return 1;
        }
        out << "threads,prefill,operations,pop_prob,min_key,max_key,"
               "seed,failed_pops,prefill_time,work_time"
            << ",cache_accesses,cache_misses" << '\n';
        out << settings.num_threads << ',' << settings.prefill_per_thread << ',' << settings.operations_per_thread
            << ',' << settings.pop_prob << ',' << settings.min_key << ',' << settings.max_key << ',' << settings.seed
            << ',' << result.num_failed_pops << ',' << std::fixed << std::setprecision(3)
            << std::chrono::duration<double>(result.prefill_time).count() << ',' << std::fixed << std::setprecision(3)
            << std::chrono::duration<double>(result.work_time).count()
#ifdef USE_PAPI
            << ',' << result.cache_accesses << ',' << result.cache_misses
#else
            << ",-1,-1"
#endif
            << std::endl;
    }
    return 0;
}
