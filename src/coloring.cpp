#include "priority_queue_factory.hpp"
#include "thread_coordination.hpp"
#ifdef PQ_MQ
#include "multiqueue/config.hpp"
#endif

#include "cxxopts.hpp"

#include <atomic>
#include <cassert>
#include <chrono>
#include <condition_variable>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

using heuristic_type = unsigned long;
using node_type = unsigned long;

using PriorityQueue = util::priority_queue_type<heuristic_type, node_type>;

struct Graph {
    std::vector<std::size_t> nodes;
    std::vector<node_type> edges;
};

struct Settings {
    std::filesystem::path graph_file;
    std::filesystem::path output;
    unsigned int num_threads = 4;
    std::uint64_t seed = 1;
    node_type starting_node = 0;
};

Settings settings;
Graph graph;

struct ThreadData {
    typename PriorityQueue::handle_type pq_handle;

    explicit ThreadData(typename PriorityQueue::handle_type h) : pq_handle(std::move(h)){};

#ifdef COUNT_STATS
    std::size_t pushed_nodes{0};
    std::size_t ignored_nodes{0};
    std::size_t extracted_nodes{0};
    std::size_t failed_extracts{0};
    std::size_t processed_nodes{0};
    std::size_t idle{0};

    inline void pushed_node(unsigned int id) {
        ++stats[id].pushed_nodes;
    }
    inline void ignored_node(unsigned int id) {
        ++stats[id].ignored_nodes;
    }
    inline void extracted_node(unsigned int id) {
        ++stats[id].extracted_nodes;
    }
    inline void extract_failed(unsigned int id) {
        ++stats[id].failed_extracts;
    }
    inline void processed_node(unsigned int id) {
        ++stats[id].processed_nodes;
    }
    inline void started_idling(unsigned int id) {
        ++stats[id].idle;
    }
#else
    inline void pushed_node(unsigned int) {
    }
    inline void ignored_node(unsigned int) {
    }
    inline void extracted_node(unsigned int) {
    }
    inline void extract_failed(unsigned int) {
    }
    inline void processed_node(unsigned int) {
    }
    inline void started_idling(unsigned int) {
    }
#endif
};

class Benchmark {
   public:
    struct Data {
        std::vector<std::atomic_int> colors;
        thread_coordination::duration_type time{};
        std::atomic_uint num_failed_pops{0};
#ifdef COUNT_STATS
        std::atomic_size_t num_locking_failed{0};
        std::atomic_size_t num_resets{0};
        std::atomic_size_t use_counts{0};
        std::atomic_size_t pushed_nodes{0};
        std::atomic_size_t ignored_nodes{0};
        std::atomic_size_t extracted_nodes{0};
        std::atomic_size_t failed_extracts{0};
        std::atomic_size_t processed_nodes{0};
        std::atomic_size_t idle{0};
#endif
        std::mutex handle_mutex;
        std::mutex idle_mutex;
        std::condition_variable idle_cv;

        explicit Data(std::size_t num_nodes) : colors(num_nodes) {
        }
    };

   private:
    static void process(ThreadData& thread_data, PriorityQueue::value_type const& retval) {
      
    }

    static bool busy_wait(ThreadData& thread_data, Data& data, std::atomic_int& num_working,
                          PriorityQueue::value_type& retval) {
        auto working = num_working.fetch_sub(1, std::memory_order_relaxed) - 1;
        while (working > 0) {
            if (thread_data.pq_handle.try_pop(retval)) {
                num_working.fetch_add(1, std::memory_order_relaxed);
                return true;
            }
            working = num_working.load(std::memory_order_relaxed);
        }
        return thread_data.pq_handle.exhaustive_check(retval);
    }

    static bool idle_wait(Data& data, std::atomic_int& num_working, int& num_busy) {
        auto l = std::unique_lock(data.idle_mutex);
        if (--num_busy == 0) {
            l.unlock();
            data.idle_cv.notify_all();
        }
        data.idle_cv.wait(l, [&]() { return num_working.load(std::memory_order_relaxed) > 0 || num_busy == 0; });
        if (num_busy == 0) {
            return false;
        }
        ++num_busy;
        l.unlock();
        num_working.fetch_add(1, std::memory_order_relaxed);
        return true;
    }

    static void work_loop(ThreadData& thread_data, Data& data, std::atomic_int& num_working, int& num_busy) {
        while (true) {
            while (true) {
                PriorityQueue::value_type retval;
                if (thread_data.pq_handle.try_pop(retval)) {
                    break;
                }
                if (busy_wait(thread_data, data, num_working, retval)) {
                    break;
                }
                if (!idle_wait(data, num_working, num_busy)) {
                    return;
                }
            }
            process(retval);
        }
    }

   public:
    static void run(thread_coordination::Context ctx, Data& data, PriorityQueue& pq) {
        std::unique_lock l{data.handle_mutex};
        ThreadData thread_data{pq.get_handle()};
        l.unlock();

        std::atomic_int num_working{0};
        ctx.execute_synchronized_timed(data.time, work_loop, thread_data, data, num_working);
        if (ctx.get_id() == 0) {
            std::clog << "done\n";
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
        ("j,threads", "The number of threads", cxxopts::value<unsigned int>(settings.num_threads), "NUMBER")
        ("p,prefill", "The prefill per thread", cxxopts::value<std::size_t>(settings.prefill_per_thread), "NUMBER")
        ("n,ops", "The number of operations per thread", cxxopts::value<std::size_t>(settings.operations_per_thread), "NUMBER")
        ("d,pop-prob", "Specify the probability of pops in [0,1]", cxxopts::value<double>(settings.pop_prob), "NUMBER")
        ("l,min", "Specify the min key", cxxopts::value<key_type>(settings.min_key), "NUMBER")
        ("m,max", "Specify the max key", cxxopts::value<key_type>(settings.max_key), "NUMBER")
        ("s,seed", "Specify the initial seed", cxxopts::value<std::uint32_t>(settings.seed), "NUMBER")
#ifdef PQ_MQ
        ("c,factor", "The factor for queues", cxxopts::value<unsigned int>(mq_config.c), "NUMBER")
        ("k,stickiness", "The stickiness period", cxxopts::value<unsigned int>(mq_config.stickiness), "NUMBER")
#endif
        ("x,no-work", "Don't perform the actual benchmark", cxxopts::value<bool>(settings.no_work))
#ifdef USE_PAPI
        ("r,pc", "Use performance counters to count L1 data cache misses", cxxopts::value<bool>(settings.use_pc))
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
    if (settings.use_pc) {
        if (int ret = PAPI_library_init(PAPI_VER_CURRENT); ret != PAPI_VER_CURRENT) {
            std::cerr << "Error initializing PAPI\n";
            return 1;
        }
        if (int ret = PAPI_thread_init((unsigned long (*)())(pthread_self)); ret != PAPI_OK) {
            std::cerr << "Error initializing threads for PAPI\n";
            return 1;
        }
        if (int ret = PAPI_query_named_event("perf::L1-DCACHE-LOAD-MISSES"); ret != PAPI_OK) {
            std::cerr << "Cannot measure L1 data cache misses\n";
            return 1;
        }
    }
#endif

#ifdef PQ_MQ
    mq_config.seed = settings.seed;
    auto pq = PriorityQueue(settings.num_threads, mq_config);
#else
    auto pq = PriorityQueue(settings.num_threads);
#endif

    std::clog << "Data structure: ";
    util::describe::describe(std::clog, pq);
    std::clog << '\n';

    Benchmark::Data benchmark_data(settings.num_threads * settings.operations_per_thread);

    thread_coordination::TaskHandle<Benchmark> task_handle(settings.num_threads, std::ref(benchmark_data),
                                                           std::ref(pq));

    task_handle.join();

    std::clog << "Prefill time (s): " << std::setprecision(3)
              << std::chrono::duration<double>(benchmark_data.prefill_time).count() << '\n';
    std::clog << "Work time (s): " << std::setprecision(3)
              << std::chrono::duration<double>(benchmark_data.work_time).count() << '\n';
    std::clog << "Failed pops: " << benchmark_data.num_failed_pops << '\n';
#ifdef USE_PAPI
    if (settings.use_pc) {
        std::clog << "L1 data cache misses: " << benchmark_data.cache_misses << '\n';
    }
#endif
#ifdef COUNT_STATS
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
    if (settings.use_pc) {
        std::cout << ',' << benchmark_data.cache_misses;
    } else {
        std::cout << ",n/a";
    }
#else
    std::cout << ",n/a";
#endif
#ifdef COUNT_STATS
    std::cout << ',' << benchmark_data.num_resets << ',' << benchmark_data.use_counts;
#else
    std::cout << ",n/a,n/a";
#endif
    std::cout << std::endl;
    return 0;
}
