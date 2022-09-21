#include "utils/priority_queue_factory.hpp"
#include "utils/thread_coordination.hpp"
#include "utils/threading.hpp"

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

using PQFactory = util::PriorityQueueFactory<key_type, value_type, true>;
using PriorityQueue = typename PQFactory::type;

using Handle = typename PriorityQueue::handle_type;

static constexpr key_type PopOp = static_cast<value_type>(std::numeric_limits<std::uint32_t>::max());

struct Settings {
    std::size_t prefill_size = 1'000'000;
    std::size_t num_operations = 100'000'000;
    unsigned int num_threads = 4;
    std::uint64_t seed = 1;
    util::PriorityQueueConfig pq_config;
    key_type min_key = 1;
    key_type max_key = PopOp >> 2;
    double pop_prob = 0.5;
};

struct alignas(2 * L1_CACHE_LINESIZE) ThreadData {
    unsigned int num_ops = 0;
    unsigned int num_failed_pops = 0;
#ifdef USE_PAPI
    long long cache_accesses;
    long long cache_misses;
#endif
};

static Settings settings;
static ThreadData* thread_data;
static key_type** prefill_keys;
static key_type* keys;
static std::chrono::steady_clock::duration prefill_time;
static std::chrono::steady_clock::duration work_time;

bool check_operations() {
    auto pop_ops = static_cast<std::size_t>(std::count(keys, keys + settings.num_operations, PopOp));
    return settings.num_threads * (settings.prefill_size / settings.num_threads) + settings.num_operations - pop_ops >=
        pop_ops;
}

struct Task {
    template <typename Generator>
    static void generate_prefill_keys(thread_coordination::Context& ctx, Generator& g) {
        auto dist = std::uniform_int_distribution<key_type>(settings.min_key, settings.max_key);
        std::size_t prefill_by_thread = settings.prefill_size / settings.num_threads;
        prefill_keys[ctx.get_id()] = new key_type[prefill_by_thread];
        std::generate(prefill_keys[ctx.get_id()], prefill_keys[ctx.get_id()] + prefill_by_thread,
                      [&]() { return dist(g); });
    }

    template <typename Generator>
    static void generate_keys(thread_coordination::Context& ctx, Generator& g) {
        auto key_dist = std::uniform_int_distribution<key_type>(settings.min_key, settings.max_key);
        auto pop_dist = std::uniform_real_distribution<>();
        ctx.execute_synchronized_blockwise(keys, keys + settings.num_operations, [&](key_type* begin, key_type* end) {
            std::generate(begin, end, [&]() { return pop_dist(g) < settings.pop_prob ? PopOp : key_dist(g); });
        });
    }

    static void prefill(thread_coordination::Context& ctx, Handle& handle) {
        ctx.execute_synchronized_timed(prefill_time, [id = ctx.get_id(), &handle]() {
            std::size_t prefill_by_thread = settings.prefill_size / settings.num_threads;
            for (std::size_t i = 0; i < prefill_by_thread; ++i) {
                handle.push({prefill_keys[id][i], prefill_keys[id][i]});
            }
        });
    }

    static void work(thread_coordination::Context& ctx, Handle& handle) {
#ifdef USE_PAPI
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
#endif
        unsigned int num_ops = 0;
        unsigned int num_failed_pops = 0;
        ctx.execute_synchronized_blockwise_timed(
            work_time, keys, keys + settings.num_operations,
            [&handle, &num_ops, &num_failed_pops](key_type* begin, const key_type* end) {
                for (auto* it = begin; it != end; ++it) {
#ifndef NULL_WORK
                    if (*it == PopOp) {
                        PriorityQueue::value_type retval;
                        // Let retval escape
                        asm volatile("" ::"g"(&retval));
                        while (!handle.try_pop(retval)) {
                            ++num_failed_pops;
                        }
                        // "Use" memory to force write to retval
                        asm volatile("" ::: "memory");
                    } else {
                        handle.push({*it, *it});
                    }
#else
                    key_type key;
                    // Let retval escape
                    asm volatile("" ::"g"(&key));
                    key = *it;
                    // "Use" memory to force write to retval
                    asm volatile("" ::: "memory");
#endif
                }
                num_ops += static_cast<unsigned int>(end - begin);
            });
#ifdef USE_PAPI
        {
            long long values[2];
            if (int ret = PAPI_stop(event_set, values); ret != PAPI_OK) {
                ctx.write(std::cerr) << "Failed to stop counters\n";
                std::abort();
            }
            thread_data[ctx.get_id()].cache_accesses = values[0];
            thread_data[ctx.get_id()].cache_misses = values[1];
        }
#endif
        thread_data[ctx.get_id()].num_ops = num_ops;
        thread_data[ctx.get_id()].num_failed_pops = num_failed_pops;
    }

    static void run(thread_coordination::Context ctx, std::uint64_t* seeds, PriorityQueue& pq) {
        if (ctx.get_id() == 0) {
            std::clog << "Generating keys..." << std::flush;
        }

        std::default_random_engine rng(seeds[ctx.get_id()]);
        generate_prefill_keys(ctx, rng);
        generate_keys(ctx, rng);

        if (ctx.get_id() == 0) {
            std::clog << "done" << std::endl;
        }

        Handle handle = ctx.execute_exclusive([&pq]() { return pq.get_handle(); });

        if (ctx.get_id() == 0) {
            std::clog << "Prefilling..." << std::flush;
        }

        prefill(ctx, handle);

        if (ctx.get_id() == 0) {
            std::clog << "done" << std::endl;
            if (!check_operations()) {
                std::cerr << "Error: More elements will be popped than pushed, "
                             "increase "
                             "prefill or decrease pop probability\n";
                std::exit(1);
            }

            std::clog << "Working..." << std::flush;
        }

        work(ctx, handle);

        if (ctx.get_id() == 0) {
            std::clog << "done\n" << std::endl;
        }
    }

    static threading::thread_config get_config(unsigned int id) {
        threading::thread_config cfg;
        cfg.cpu_set.reset();
        cfg.cpu_set.set(id);
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

    std::filesystem::path out_file;

    cxxopts::Options options(argv[0]);
    options.add_options()
        // clang-format off
        ("h,help", "Print this help")
        ("p,prefill", "Specify the number of elements to prefill the queue with "
         "(default: 1'000'000)", cxxopts::value<size_t>(), "NUMBER")
        ("n,ops", "Specify the number of operations"
         "(default: 100'000'000)", cxxopts::value<std::size_t>(), "NUMBER")
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
        ("k,stickiness", "The stickiness (only has effect for some multiqueue variants)"
         "(default: 8)", cxxopts::value<unsigned int>(), "NUMBER")
        ("o,output", "Output data in csv (comma-separated)", cxxopts::value<std::filesystem::path>(out_file), "PATH")
        // clang-format on
        ;

    try {
        auto result = options.parse(argc, argv);
        if (result.count("help") > 0) {
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
        if (result.count("stickiness") > 0) {
            settings.pq_config.stickiness = result["stickiness"].as<unsigned int>();
        }
    } catch (cxxopts::OptionParseException const& e) {
        std::cerr << "Error parsing arguments: " << e.what() << '\n';
        std::cerr << options.help() << std::endl;
        return 1;
    }

    std::cout << "Settings\n"
              << "Prefill: " << settings.prefill_size << '\n'
              << "Operations: " << settings.num_operations << '\n'
              << "Threads: " << settings.num_threads << '\n'
              << "Pop probability: " << std::fixed << settings.pop_prob << '\n'
              << "Min key: " << settings.min_key << '\n'
              << "Max key: " << settings.max_key << '\n'
              << "Seed: " << settings.seed << '\n';
    std::cout << '\n';

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

    std::default_random_engine rng(settings.seed);

    settings.pq_config.seed = rng();
    auto pq = PQFactory::create(settings.num_threads, settings.pq_config);

    std::cout << "Priority queue\nName: ";
    PQFactory::describe(std::cout, pq);
    std::cout << '\n';

    thread_data = new ThreadData[settings.num_threads];
    prefill_keys = new key_type*[settings.num_threads];
    keys = new key_type[settings.num_operations];
    auto thread_seeds = std::make_unique<std::uint64_t[]>(settings.num_threads);
    for (std::size_t i = 0; i < settings.num_threads; ++i) {
        thread_seeds[i] = rng();
    }
    thread_coordination::run_task<Task>(settings.num_threads, thread_seeds.get(), std::ref(pq));
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
