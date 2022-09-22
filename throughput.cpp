#include "utils/priority_queue_factory.hpp"
#include "utils/thread_coordination.hpp"
#include "utils/threading.hpp"
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

using PriorityQueue = util::priority_queue_type<key_type, value_type, true>;
using Handle = typename PriorityQueue::handle_type;

static constexpr key_type PopOp = static_cast<value_type>(std::numeric_limits<std::uint32_t>::max());
static constexpr std::size_t DefaultPrefillPerThread = 1 << 20;
static constexpr std::size_t DefaultOperationsPerThread = 1 << 24;
static constexpr int DefaultNumThreads = 4;
static constexpr double DefaultPopProbability = 0.5;

struct Settings {
    std::size_t prefill_per_thread = DefaultPrefillPerThread;
    std::size_t operations_per_thread = DefaultOperationsPerThread;
    int num_threads = DefaultNumThreads;
    double pop_prob = DefaultPopProbability;
    std::uint32_t seed = 1;
    key_type min_key = 1;
    key_type max_key = PopOp >> 2;
};

struct Stats {
    thread_coordination::duration_type prefill_time;
    thread_coordination::duration_type work_time;
    std::atomic_size_t min_num_ops = std::numeric_limits<std::size_t>::max();
    std::atomic_size_t max_num_ops = 0;
    std::atomic_size_t num_failed_pops = 0;
#ifdef USE_PAPI
    std::atomic_size_t cache_accesses = 0;
    std::atomic_size_t cache_misses = 0;
#endif
};

bool check_operations(Settings const& settings, std::vector<key_type> const& keys) {
    std::size_t num_elements = static_cast<std::size_t>(settings.num_threads) * settings.prefill_per_thread;
    for (auto k : keys) {
        if (k == PopOp) {
            if (num_elements == 0) {
                return false;
            }
            --num_elements;
        } else {
            ++num_elements;
        }
    }
    return true;
}

struct Task {
    struct ThreadData {
        std::size_t num_ops = 0;
        std::size_t num_failed_pops = 0;
#ifdef USE_PAPI
        long long cache_accesses = 0;
        long long cache_misses = 0;
#endif
    };

    static std::mutex m;
    static std::condition_variable cv;
    static bool check_done;
    static bool check_successful;

    template <typename Generator>
    static std::vector<key_type> generate_prefill(Settings const& settings, Generator& g) {
        auto key_dist = std::uniform_int_distribution<key_type>(settings.min_key, settings.max_key);
        std::vector<key_type> keys(settings.prefill_per_thread);
        std::generate(keys.begin(), keys.end(), [&]() { return key_dist(g); });
        return keys;
    }

    template <typename Generator>
    static void fill_keys(thread_coordination::Context const& ctx, Settings const& settings,
                          std::vector<key_type>& keys, Generator& g) {
        auto key_dist = std::uniform_int_distribution<key_type>(settings.min_key, settings.max_key);
        auto pop_dist = std::bernoulli_distribution(settings.pop_prob);
        ctx.execute_synchronized_blockwise(keys.begin(), keys.end(), [&](auto block_begin, auto block_end) {
            std::generate(block_begin, block_end, [&]() { return pop_dist(g) ? PopOp : key_dist(g); });
        });
    }

    static void prefill(thread_coordination::Context const& ctx, Handle& handle, std::vector<key_type> const& keys,
                        Stats& stats) {
        ctx.execute_synchronized_timed(stats.prefill_time, [&handle, &keys]() {
            for (auto k : keys) {
                handle.push({k, k});
            }
        });
    }

#ifdef USE_PAPI
    static void start_performance_counter() {
        if (int ret = PAPI_register_thread(); ret != PAPI_OK) {
            ctx.write(std::cerr) << "Failed to register thread to PAPI\n";
            std::abort();
        }
        int event_set = PAPI_NULL;
        if (int ret = PAPI_create_eventset(&event_set); ret != PAPI_OK) {
            ctx.write(std::cerr) << "Failed to init event set\n";
            std::abort();
        }
        if (int ret = PAPI_add_event(event_set, PAPI_L1_DCA); ret != PAPI_OK) {
            ctx.write(std::cerr) << "Failed to count cache accesses\n";
            std::abort();
        }
        if (int ret = PAPI_add_named_event(event_set, "perf::L1-DCACHE-LOAD-MISSES"); ret != PAPI_OK) {
            ctx.write(std::cerr) << "Failed to count cache misses\n";
            std::abort();
        }
        if (int ret = PAPI_start(event_set); ret != PAPI_OK) {
            ctx.write(std::cerr) << "Failed to start counters\n";
            std::abort();
        }
    }
#endif

    static void work(thread_coordination::Context const& ctx, Handle& handle, std::vector<key_type> const& keys,
                     Stats& stats) {
        ThreadData data;
#ifdef USE_PAPI
        start_performance_counter();
#endif
        ctx.execute_synchronized_blockwise_timed(stats.work_time, keys.begin(), keys.end(),
                                                 [&handle, &data](auto block_begin, auto block_end) {
                                                     for (auto it = block_begin; it != block_end; ++it) {
                                                         if (*it == PopOp) {
                                                             PriorityQueue::value_type retval;
                                                             // Let retval escape
                                                             asm volatile("" ::"g"(&retval));
                                                             while (!handle.try_pop(retval)) {
                                                                 ++data.num_failed_pops;
                                                             }
                                                             // "Use" memory to force write to retval
                                                             asm volatile("" ::: "memory");
                                                         } else {
                                                             handle.push({*it, *it});
                                                         }
                                                     }
                                                     data.num_ops += static_cast<std::size_t>(block_end - block_begin);
                                                 });
#ifdef USE_PAPI
        {
            long long values[2];
            if (int ret = PAPI_stop(event_set, values); ret != PAPI_OK) {
                ctx.write(std::cerr) << "Failed to stop counters\n";
                std::abort();
            }
            data.cache_accesses = values[0];
            data.cache_misses = values[1];
        }
#endif
    }

    static void run(thread_coordination::Context ctx, Settings const& settings, std::vector<key_type>& keys,
                    Stats& stats, PriorityQueue& pq) {
        if (ctx.get_id() == 0) {
            std::clog << "Generating keys..." << std::flush;
            check_done = false;
        }
        std::seed_seq seed{settings.seed, static_cast<std::uint32_t>(ctx.get_id())};
        std::default_random_engine rng(seed);
        std::vector<key_type> prefill_keys = generate_prefill(settings, rng);
        fill_keys(ctx, settings, keys, rng);

        if (ctx.get_id() == 0) {
            std::clog << "done" << std::endl;
            check_successful = check_operations(settings, keys);
            {
                std::unique_lock l(m);
                check_done = true;
            }
            cv.notify_all();
        } else {
            std::unique_lock l(m);
            cv.wait(l, [] { return check_done; });
        }

        Handle handle = ctx.execute_exclusive([&pq]() { return pq.get_handle(); });

        if (ctx.get_id() == 0) {
            std::clog << "Prefilling..." << std::flush;
        }

        prefill(ctx, handle, prefill_keys, stats);

        if (ctx.get_id() == 0) {
            std::clog << "done\nWorking..." << std::flush;
        }

        work(ctx, handle, keys, stats);

        if (ctx.get_id() == 0) {
            std::clog << "done\n" << std::endl;
        }
    }

    static threading::thread_config get_config(int id) {
        threading::thread_config cfg;
        cfg.cpu_set.set(static_cast<std::size_t>(id));
        return cfg;
    }
};

int main(int argc, char* argv[]) {
    std::cout << "Build config\n";
#ifndef NDEBUG
    std::cout << "DEBUG: enabled\n";
#else
    std::cout << "DEBUG: disabled\n";
#endif
#ifdef USE_PAPI
    std::cout << "PAPI: enabled\n";
#else
    std::cout << "PAPI: disabled\n";
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
        ("p,prefill", "The prefill per thread", cxxopts::value<std::size_t>(settings.prefill_per_thread), "NUMBER")
        ("n,ops", "The number of operations per thread", cxxopts::value<std::size_t>(settings.operations_per_thread), "NUMBER")
        ("j,threads", "The number of threads", cxxopts::value<int>(settings.num_threads), "NUMBER")
        ("d,pop-prob", "Specify the probability of pops in [0,1]", cxxopts::value<double>(settings.pop_prob), "NUMBER")
        ("l,min", "Specify the min key", cxxopts::value<key_type>(settings.min_key), "NUMBER")
        ("m,max", "Specify the max key", cxxopts::value<key_type>(settings.max_key), "NUMBER")
        ("s,seed", "Specify the initial seed", cxxopts::value<std::uint32_t>(settings.seed), "NUMBER")
#ifdef PQ_MQ
        ("c,factor", "The factor for queues", cxxopts::value<int>(mq_config.c), "NUMBER")
        ("k,stickiness", "The stickiness period", cxxopts::value<int>(mq_config.stickiness), "NUMBER")
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

    std::cout << "Settings\n"
              << "Prefill: " << settings.prefill_per_thread << '\n'
              << "Operations: " << settings.operations_per_thread << '\n'
              << "Threads: " << settings.num_threads << '\n'
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
    if (int ret = PAPI_thread_init((unsigned long (*)(void))(pthread_self)); ret != PAPI_OK) {
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

    std::vector<key_type> keys(settings.operations_per_thread * static_cast<std::size_t>(settings.num_threads));

    auto thread_seeds = std::make_unique<std::uint64_t[]>(settings.num_threads);
    for (std::size_t i = 0; i < settings.num_threads; ++i) {
        thread_seeds[i] = rng();
    }
    thread_coordination::run_task<Task>(settings.num_threads, thread_seeds.get(), std::ref(pq));
    if (Task::check_successful) {
        std::cerr << "Error: More elements will be popped than pushed, "
                     "increase "
                     "prefill or decrease pop probability\n";
    }

    unsigned int total_failed_pops =
        std::accumulate(thread_data, thread_data + settings.num_threads, 0U,
                        [](unsigned int sum, auto const& d) { return sum + d.num_failed_pops; });
    auto [min_thread, max_thread] =
        std::minmax_element(thread_data, thread_data + settings.num_threads,
                            [](auto const& lhs, auto const& rhs) { return lhs.num_ops < rhs.num_ops; });
#ifdef USE_PAPI
    auto cache_accesses = std::accumulate(thread_data, thread_data + settings.num_threads, .0,
                                          [](auto sum, auto const& s) { return sum + s.cache_accesses; });
    auto cache_misses = std::accumulate(thread_data, thread_data + settings.num_threads, .0,
                                        [](auto sum, auto const& s) { return sum + s.cache_misses; });
#endif
    std::cout << "Least/Most operations by any thread: " << min_thread->num_ops << '/' << max_thread->num_ops << '\n';
    std::cout << "Failed pops: " << total_failed_pops << '\n';
    std::cout << "Prefill time (s): " << std::setprecision(3) << std::chrono::duration<double>(prefill_time).count()
              << '\n';
    std::cout << "Work time (s): " << std::setprecision(3) << std::chrono::duration<double>(work_time).count() << '\n';
#ifdef USE_PAPI
    std::cout << "Cache accesses/misses: " << cache_accesses << '/' << cache_misses << '\n';
#endif

    if (!out_file.empty()) {
        auto out = std::ofstream{out_file};
        if (!out.is_open()) {
            std::cerr << "Could not open file to write out results\n";
            return 1;
        }
        out << "prefill,ops,threads,pop_prob,min_key,max_key,"
               "seed,min_thread_ops,max_thread_ops,failed_pops,prefill_time,"
               "work_time"
#ifdef USE_PAPI
            << ",cache_accesses,cache_misses"
#endif
            << '\n';
        out << settings.prefill_size << ',' << settings.num_operations << ',' << settings.num_threads << ','
            << settings.pop_prob << ',' << settings.min_key << ',' << settings.max_key << ',' << settings.seed << ','
            << min_thread->num_ops << ',' << max_thread->num_ops << ',' << total_failed_pops << ',' << std::fixed
            << std::setprecision(2) << std::chrono::duration<double>(prefill_time).count() << ',' << std::fixed
            << std::setprecision(2) << std::chrono::duration<double>(work_time).count()
#ifdef USE_PAPI
            << ',' << cache_accesses << ',' << cache_misses
#endif
            << std::endl;
    }

    for (unsigned int i = 0; i < settings.num_threads; ++i) {
        delete[] prefill_keys[i];
    }
    delete[] prefill_keys;
    delete[] keys;
    delete[] thread_data;
    return 0;
}
