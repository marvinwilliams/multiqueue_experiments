#ifndef L1_CACHE_LINESIZE
#define L1_CACHE_LINESIZE 64
#endif

#include "utils/priority_queue_factory.hpp"
#include "utils/thread_coordination.hpp"
#include "utils/threading.hpp"
#include "utils/xoroshiro256starstar.hpp"

#include "external/cxxopts.hpp"

#ifdef USE_PAPI
#include <papi.h>
#include <pthread.h>
#endif

#include <atomic>
#include <cassert>
#include <chrono>
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

using PriorityQueue =
    typename util::PriorityQueueFactory<key_type, value_type, true>::type;
static constexpr key_type PopOp =
    static_cast<value_type>(std::numeric_limits<std::uint32_t>::max());

struct Settings {
    std::size_t prefill_size = 1'000'000;
    std::size_t num_operations = 10'000'000;
    unsigned int num_threads = 4;
    std::uint64_t seed = 1;
    util::PriorityQueueConfig pq_config;
    key_type min_key = 1;
    key_type max_key = PopOp >> 2;
    double pop_prob = 0.5;
};

struct alignas(L1_CACHE_LINESIZE) ThreadData {
    xoroshiro256starstar rng;
    unsigned int num_ops = 0;
    unsigned int num_failed_pops = 0;
#ifdef USE_PAPI
    long long l1_stats[2];
#endif
};

static Settings settings;
alignas(L1_CACHE_LINESIZE) static ThreadData* thread_data;
alignas(L1_CACHE_LINESIZE) static key_type* prefill_keys;
alignas(L1_CACHE_LINESIZE) static key_type* keys;

void generate_prefill_keys(thread_coordination::Context& ctx) {
    auto dist = std::uniform_int_distribution<key_type>(settings.min_key,
                                                        settings.max_key);
    ctx.execute_synchronized_blockwise(
        prefill_keys, prefill_keys + settings.prefill_size,
        [id = ctx.get_id(), &dist](key_type* begin, key_type* end) {
            std::generate(begin, end,
                          [id, &dist]() { return dist(thread_data[id].rng); });
        });
}

void generate_keys(thread_coordination::Context& ctx) {
    auto dist = std::uniform_int_distribution<key_type>(settings.min_key,
                                                        settings.max_key);
    auto pop_dist = std::uniform_real_distribution<>();
    ctx.execute_synchronized_blockwise(
        keys, keys + settings.num_threads * settings.num_operations,
        [id = ctx.get_id(), &dist, &pop_dist](key_type* begin, key_type* end) {
            std::generate(begin, end, [id, &dist, &pop_dist]() {
                return pop_dist(thread_data[id].rng) < settings.pop_prob
                           ? PopOp
                           : dist(thread_data[id].rng);
            });
        });
}

bool check_operations() {
    auto pop_ops = static_cast<std::size_t>(std::count(
        keys, keys + settings.num_threads * settings.num_operations, PopOp));
    return settings.prefill_size +
               settings.num_threads * settings.num_operations - pop_ops >=
           pop_ops;
}

void prefill(thread_coordination::Context& ctx, PriorityQueue::Handle& handle,
             std::chrono::steady_clock::duration& duration) {
    ctx.execute_synchronized_blockwise_timed(
        duration, prefill_keys, prefill_keys + settings.prefill_size,
        [&handle, id = ctx.get_id()](key_type* begin, key_type* end) {
            for (auto it = begin; it != end; ++it) {
                assert(*it < PopOp);
                handle.push({*it, *it});
            };
        });
}

void work(thread_coordination::Context& ctx, PriorityQueue::Handle& handle,
          std::chrono::steady_clock::duration& duration) {
#ifdef USE_PAPI
    int event_set = PAPI_NULL;
    if (int ret = PAPI_create_eventset(&event_set); ret != PAPI_OK) {
        ctx.write(std::cerr) << "Failed to init event set\n";
        std::abort();
    }
    if (int ret = PAPI_add_event(event_set, PAPI_L1_DCA); ret != PAPI_OK) {
        ctx.write(std::cerr) << "Failed to add PAPI_L1_DCA to event set";
        std::abort();
    }
    if (int ret = PAPI_add_event(event_set, PAPI_BR_MSP); ret != PAPI_OK) {
        ctx.write(std::cerr) << "Failed to add PAPI_BR_MSP to event set\n";
        std::abort();
    }
    if (int ret = PAPI_start(event_set); ret != PAPI_OK) {
        ctx.write(std::cerr) << "Failed to start counters\n";
        std::abort();
    }
#endif
    ctx.execute_synchronized_blockwise_timed(
        duration, keys, keys + settings.num_threads * settings.num_operations,
        [&, id = ctx.get_id()](key_type* begin, key_type* end) {
            for (auto it = begin; it != end; ++it) {
                if (*it == PopOp) {
                    PriorityQueue::value_type retval;
                    // Let retval escape
                    asm volatile("" ::"g"(&retval));
                    while (!handle.try_pop(retval)) {
                        ++thread_data[id].num_failed_pops;
                    }
                    // "Use" memory to force write to retval
                    asm volatile("" ::: "memory");
                } else {
                    handle.push({*it, *it});
                }
            }
            thread_data[id].num_ops += static_cast<unsigned int>(end - begin);
        });
#ifdef USE_PAPI
    if (int ret = PAPI_stop(event_set, thread_data[ctx.get_id()].l1_stats);
        ret != PAPI_OK) {
        ctx.write(std::cerr) << "Failed to stop counters\n";
        std::abort();
    }
#endif
}

struct GenerateTask {
    static void run(thread_coordination::Context ctx) {
        generate_prefill_keys(ctx);
        generate_keys(ctx);
    }

    static threading::thread_config get_config(unsigned int id) {
        threading::thread_config cfg;
        cfg.cpu_set.reset();
        cfg.cpu_set.set(id);
        return cfg;
    }
};

struct WorkTask {
    static void run(thread_coordination::Context ctx, PriorityQueue& pq,
                    std::chrono::steady_clock::duration& prefill_duration,
                    std::chrono::steady_clock::duration& work_duration) {
#ifdef USE_PAPI
        if (int ret = PAPI_register_thread(); ret != PAPI_OK) {
            ctx.write(std::cerr) << "Failed to register thread to PAPI\n";
            std::abort();
        }
#endif
        PriorityQueue::Handle handle = pq.get_handle();

        prefill(ctx, handle, prefill_duration);

        if (ctx.get_id() == 0) {
            std::clog << "done\nWorking..." << std::flush;
        }
        work(ctx, handle, work_duration);
    }

    static threading::thread_config get_config(unsigned int id) {
        threading::thread_config cfg;
        cfg.cpu_set.reset();
        cfg.cpu_set.set(id);
        return cfg;
    }
};

int main(int argc, char* argv[]) {
#ifndef NDEBUG
    std::clog << "DEBUG build\n";
#endif
    std::clog << "L1 cache linesize (byte): " << L1_CACHE_LINESIZE << "\n\t";
    std::clog << '\n';

    std::clog << "Command line: ";
    std::copy(argv, argv + argc,
              std::ostream_iterator<char const*>(std::clog, " "));
    std::clog << "\n\n";

    cxxopts::Options options(argv[0]);
    options.add_options()
        // clang-format off
        ("h,help", "Print this help")
        ("p,prefill", "Specify the number of elements to prefill the queue with "
         "(default: 1'000'000)", cxxopts::value<size_t>(), "NUMBER")
        ("n,ops", "Specify the number of operations per thread"
         "(default: 10'000'000)", cxxopts::value<std::size_t>(), "NUMBER")
        ("j,threads", "Specify the number of threads "
         "(default: 4)", cxxopts::value<unsigned int>(), "NUMBER")
        ("d,pop-prob", "Specify the probability of pops"
         "(default: 0.5)", cxxopts::value<double>(), "NUMBER")
        ("m,max", "Specify the max key "
         "(default: MAX)", cxxopts::value<key_type>(), "NUMBER")
        ("l,min", "Specify the min key "
         "(default: 0)", cxxopts::value<key_type>(), "NUMBER")
        ("s,seed", "Specify the initial seed"
         "(default: 0)", cxxopts::value<std::uint32_t>(), "NUMBER")
        ("c,factor", "The number of queues when using multiqueue or multififo"
         "(default: 4)", cxxopts::value<std::size_t>(), "NUMBER")
#ifdef MQ_HAS_STICKINESS
        ("k,stickiness", "The stickiness"
         "(default: 8)", cxxopts::value<unsigned int>(), "NUMBER")
#endif
        // clang-format on
        ;

    try {
        auto result = options.parse(argc, argv);
        if (result.count("help")) {
            std::cerr << options.help() << std::endl;
            return 0;
        }
        if (result.count("prefill") > 0) {
            settings.prefill_size = result["prefill"].as<size_t>();
        }
        if (result.count("ops") > 0) {
            settings.num_operations = result["ops"].as<std::size_t>();
        }
        if (result.count("threads") > 0) {
            settings.num_threads = result["threads"].as<unsigned int>();
        }
        if (result.count("pop-prob") > 0) {
            settings.pop_prob = result["pop-prob"].as<double>();
        }
        if (result.count("max") > 0) {
            settings.max_key = result["max"].as<key_type>();
        }
        if (result.count("min") > 0) {
            settings.min_key = result["min"].as<key_type>();
        }
        if (result.count("seed") > 0) {
            settings.seed = result["seed"].as<std::uint32_t>();
        }
        if (result.count("factor") > 0) {
            settings.pq_config.c = result["factor"].as<std::size_t>();
        }
#ifdef MQ_HAS_STICKINESS
        if (result.count("stickiness") > 0) {
            settings.pq_config.stickiness =
                result["stickiness"].as<unsigned int>();
        }
#endif
    } catch (cxxopts::OptionParseException const& e) {
        std::cerr << "Error parsing arguments: " << e.what() << '\n';
        std::cerr << options.help() << std::endl;
        return 1;
    }

    std::clog << "Settings: \n\t"
              << "Prefill: " << settings.prefill_size << "\n\t"
              << "Operations per thread: " << settings.num_operations << "\n\t"
              << "Threads: " << settings.num_threads << "\n\t"
              << "Pop probability: " << std::fixed << settings.pop_prob
              << "\n\t"
              << "Max key: " << settings.max_key << "\n\t"
              << "Min key: " << settings.min_key << "\n\t"
              << "Seed: " << settings.seed;
    std::clog << "\n\n";

    xoroshiro256starstar rng(settings.seed);

    settings.pq_config.seed = rng();
    auto pq = util::PriorityQueueFactory<key_type, value_type, true>::create(
        settings.num_threads, settings.pq_config);

    thread_data = ::new ThreadData[settings.num_threads];
    for (std::size_t i = 0; i < settings.num_threads; ++i) {
        thread_data[i].rng.seed(rng());
    }

    prefill_keys = ::new (std::align_val_t{L1_CACHE_LINESIZE})
        key_type[settings.prefill_size];

    keys = ::new (std::align_val_t{L1_CACHE_LINESIZE})
        key_type[settings.num_threads * settings.num_operations];
    std::clog << "Generating keys..." << std::flush;
    thread_coordination::run_task<GenerateTask>(settings.num_threads);
    std::clog << "done\n";

    if (!check_operations()) {
        std::cerr << "More elements will be popped than pushed, increase "
                     "prefill or decrease pop probability\n";
        return 1;
    }

#ifdef USE_PAPI
    if (int ret = PAPI_library_init(PAPI_VER_CURRENT);
        ret != PAPI_VER_CURRENT) {
        std::cerr << "Error initializing PAPI\n";
        return 1;
    }
    if (int ret = PAPI_thread_init((unsigned long (*)(void))(pthread_self));
        ret != PAPI_OK) {
        std::cerr << "Error initializing threads for PAPI\n";
        return 1;
    }
    if (int ret = PAPI_query_event(PAPI_L1_DCA); ret != PAPI_OK) {
        std::cerr << "PAPI_L1_DCA not available\n";
        return 1;
    }
    if (int ret = PAPI_query_event(PAPI_BR_MSP); ret != PAPI_OK) {
        std::cerr << "PAPI_BR_MSP not available\n";
        return 1;
    }
#endif

    std::clog << "Prefilling... " << std::flush;
    std::chrono::steady_clock::duration prefill_duration;
    std::chrono::steady_clock::duration work_duration;
    thread_coordination::run_task<WorkTask>(settings.num_threads, std::ref(pq),
                                            std::ref(prefill_duration),
                                            std::ref(work_duration));
    std::clog << "done\n\n";

    unsigned int total_failed_pops =
        std::accumulate(thread_data, thread_data + settings.num_threads, 0u,
                        [](unsigned int sum, auto const& data) {
                            return sum + data.num_failed_pops;
                        });
    auto [min_thread, max_thread] =
        std::minmax_element(thread_data, thread_data + settings.num_threads,
                            [](auto const& lhs, auto const& rhs) {
                                return lhs.num_ops < rhs.num_ops;
                            });
#ifdef USE_PAPI
    auto avg_cache_accesses = std::accumulate(
        thread_data, thread_data + settings.num_threads, .0,
        [](double avg, auto const& data) {
            return avg + static_cast<double>(data.l1_stats[0]) /
                             static_cast<double>(settings.num_threads *
                                                 settings.num_operations);
        });
    auto avg_cache_misses = std::accumulate(
        thread_data, thread_data + settings.num_threads, .0,
        [](double avg, auto const& data) {
            return avg + static_cast<double>(data.l1_stats[1]) /
                             static_cast<double>(settings.num_threads *
                                                 settings.num_operations);
        });
#endif
    std::clog << "Most operations: " << max_thread->num_ops << '\n';
    std::clog << "Least operations: " << min_thread->num_ops << '\n';
    std::clog << "Failed pops: " << total_failed_pops << '\n';
    std::clog << "Prefill time (s): " << std::setprecision(3)
              << std::chrono::duration<double>(prefill_duration).count()
              << '\n';
    std::clog << "Work time (s): " << std::setprecision(3)
              << std::chrono::duration<double>(work_duration).count() << '\n';
    std::clog << "Throughput per thread: " << std::fixed << std::setprecision(2)
              << static_cast<double>(settings.num_threads *
                                     settings.num_operations) /
                     std::chrono::duration<double>(work_duration).count()
              << '\n';
#ifdef USE_PAPI
    std::clog << "Avg cache accesses: " << avg_cache_accesses << '\n';
    std::clog << "Avg cache misses: " << avg_cache_misses << '\n';
#endif

    std::cout << "prefill,ops_per_thread,threads,pop_prob,min_key,max_key,"
                 "seed,min_ops,max_ops,failed_pops,prefill_time,work_time,"
                 "throughput\n"
              << settings.prefill_size << ',' << settings.num_operations << ','
              << settings.num_threads << ',' << settings.pop_prob << ','
              << settings.min_key << ',' << settings.max_key << ','
              << settings.seed << ',' << min_thread->num_ops << ','
              << max_thread->num_ops << ',' << total_failed_pops << ','
              << std::fixed << std::setprecision(2)
              << std::chrono::duration<double>(prefill_duration).count() << ','
              << std::fixed << std::setprecision(2)
              << std::chrono::duration<double>(work_duration).count() << ','
              << std::fixed << std::setprecision(2)
              << static_cast<double>(settings.num_threads *
                                     settings.num_operations) /
                     std::chrono::duration<double>(work_duration).count()
              << std::endl;

    ::operator delete[](keys, std::align_val_t{L1_CACHE_LINESIZE});
    ::operator delete[](prefill_keys, std::align_val_t{L1_CACHE_LINESIZE});
    delete[] thread_data;
    return 0;
}
