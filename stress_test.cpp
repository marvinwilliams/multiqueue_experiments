#include "utils/operation_generator.hpp"
#include "utils/priority_queue_factory.hpp"
#include "utils/thread_coordination.hpp"
#include "utils/threading.hpp"
#include "utils/xoroshiro256starstar.hpp"

#include "external/cxxopts.hpp"

#include <time.h>
#include <x86intrin.h>
#include <atomic>
#include <cassert>
#include <chrono>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <new>
#include <numeric>
#include <optional>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

#if !defined THROUGHPUT_MODE && !defined QUALITY_MODE
#error Need to define either THROUGHPUT_MODE or QUALITY_MODE
#endif

#if defined THROUGHPUT_MODE
#undef QUALITY_MODE
#endif
#if defined QUALITY_MODE
#undef THROUGHPUT_MODE
#endif

#ifndef L1_CACHE_LINESIZE
#define L1_CACHE_LINESIZE 64
#endif

#ifndef PAGESIZE
#define PAGESIZE 4096
#endif

using key_type = unsigned long;
using value_type = unsigned long;

using PriorityQueue =
    typename util::PriorityQueueFactory<key_type, value_type>::type;

using steady_clock = std::chrono::steady_clock;

#ifdef QUALITY_MODE

using tick_type = std::uint64_t;

static inline tick_type get_tick() noexcept {
#ifdef USE_TSC
    // Not synchronized among sockets
    _mm_lfence();
    return __rdtsc();
    _mm_lfence();
#else
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return static_cast<tick_type>(ts.tv_sec * 1000000000 + ts.tv_nsec);
#endif
}

#endif

struct Settings {
    std::size_t prefill_size = 1'000'000;
    std::size_t num_operations = 10'000'000;
#ifdef QUALITY_MODE
    std::chrono::microseconds sleep_between_operations =
        std::chrono::microseconds::zero();
#endif
    unsigned int num_threads = 4;
    std::uint64_t seed = 1;
    util::PriorityQueueParameters pq_params;
    OperationGenerator<key_type> insert_config{
        InsertPolicy::Uniform,
        KeyDistribution::Uniform,
        1,
        static_cast<value_type>(std::numeric_limits<std::uint32_t>::max() >>
                                2),  // Some safety for sentinels
        1,
        100,
    };
};

struct OperationCount {
    std::size_t num_prefill_insertions = 0;
    std::size_t num_insertions = 0;
    std::size_t num_deletions = 0;
    std::size_t num_failed_deletions = 0;
};

OperationCount& operator+=(OperationCount& lhs,
                           OperationCount const& rhs) noexcept {
    lhs.num_prefill_insertions += rhs.num_prefill_insertions;
    lhs.num_insertions += rhs.num_insertions;
    lhs.num_deletions += rhs.num_deletions;
    lhs.num_failed_deletions += rhs.num_failed_deletions;
    return lhs;
}

OperationCount operator+(OperationCount lhs,
                         OperationCount const& rhs) noexcept {
    return lhs += rhs;
}

#ifdef QUALITY_MODE

static constexpr unsigned int bits_for_thread_id = 8;
static constexpr value_type value_mask =
    (static_cast<value_type>(1)
     << (std::numeric_limits<value_type>::digits - bits_for_thread_id)) -
    1;

static constexpr value_type to_value(unsigned int thread_id,
                                     value_type elem_id) noexcept {
    return (static_cast<value_type>(thread_id)
            << (std::numeric_limits<value_type>::digits - bits_for_thread_id)) |
           (elem_id & value_mask);
}

static constexpr unsigned int get_thread_id(value_type value) noexcept {
    return static_cast<unsigned int>(
        value >>
        (std::numeric_limits<value_type>::digits - bits_for_thread_id));
}

static constexpr value_type get_elem_id(value_type value) noexcept {
    return value & value_mask;
}

struct InsertionLogEntry {
    tick_type tick;
    key_type key;
};

struct DeletionLogEntry {
    tick_type tick;
    value_type value;
};

#endif

struct alignas(L1_CACHE_LINESIZE) ThreadData {
    xoroshiro256starstar rng;
    OperationCount op_count;
    InsertingStrategy<key_type> inserter;
#ifdef QUALITY_MODE
    std::vector<InsertionLogEntry> ins_log;
    std::vector<DeletionLogEntry> del_log;
    std::vector<tick_type> failed_del_log;
#endif
};

static Settings settings;
alignas(L1_CACHE_LINESIZE) static ThreadData* thread_data;
alignas(L1_CACHE_LINESIZE) static key_type* prefill_keys;
alignas(L1_CACHE_LINESIZE) static key_type* keys;

void generate_prefill_keys(thread_coordination::Context& ctx) {
    ctx.execute_synchronized_blockwise_timed(
        "generate_prefill_keys", prefill_keys,
        prefill_keys + settings.prefill_size,
        [id = ctx.get_id()](key_type* begin, key_type* end) {
            std::generate(begin, end, [id]() {
                return thread_data[id].inserter.get_key(thread_data[id].rng);
            });
        });
}

void generate_keys(thread_coordination::Context& ctx) {
    ctx.execute_synchronized_blockwise_timed(
        "generate_keys", keys,
        keys + settings.num_threads * settings.num_operations,
        [id = ctx.get_id()](key_type* begin, key_type* end) {
            std::generate(begin, end, [id]() {
                return thread_data[id].inserter.insert(thread_data[id].rng)
                           ? thread_data[id].inserter.get_key(
                                 thread_data[id].rng)
                           : std::numeric_limits<std::uint32_t>::max();
            });
        });
}

void prefill(thread_coordination::Context& ctx, PriorityQueue::Handle& handle) {
    ctx.execute_synchronized_blockwise_timed(
        "prefill", prefill_keys, prefill_keys + settings.prefill_size,
        [&, id = ctx.get_id()](key_type* begin, key_type* end) {
            std::for_each(begin, end, [&, id](key_type k) {
#ifdef QUALITY_MODE
                auto v = to_value(id, thread_data[id].ins_log.size());
                handle.push({k, v});
                thread_data[id].ins_log.push_back({0, k});
#else
          handle.push({k, k});
#endif
                ++thread_data[id].op_count.num_prefill_insertions;
            });
        });
}

#ifdef QUALITY_MODE

void work(thread_coordination::Context& ctx, PriorityQueue::Handle& handle) {
    ctx.execute_synchronized_blockwise_timed(
        "work", keys, keys + settings.num_threads * settings.num_operations,
        [&, id = ctx.get_id()](key_type* begin, key_type* end) {
            std::for_each(begin, end, [&, id](key_type k) {
                if (k == std::numeric_limits<std::uint32_t>::max()) {
                    PriorityQueue::value_type retval;
                    while (!handle.try_extract_top(retval)) {
                        // Not expected due to prefill
                        auto tick = get_tick();
                        thread_data[id].failed_del_log.push_back(tick);
                        ++thread_data[id].op_count.num_failed_deletions;
                    }
                    auto tick = get_tick();

                    thread_data[id].del_log.push_back({tick, retval.second});
                    ++thread_data[id].op_count.num_deletions;
                } else {
                    value_type value = to_value(
                        id, thread_data[id].op_count.num_prefill_insertions +
                                thread_data[id].op_count.num_insertions);
                    /* auto tick = get_tick(); */
                    handle.push({k, value});
                    auto tick = get_tick();
                    thread_data[id].ins_log.push_back({tick, k});
                    ++thread_data[id].op_count.num_insertions;
                }
                if (settings.sleep_between_operations !=
                    std::chrono::microseconds::zero()) {
                    auto dist = std::uniform_int_distribution<std::uint64_t>(
                        1, static_cast<std::uint64_t>(
                               settings.sleep_between_operations.count()));
                    std::this_thread::sleep_for(
                        std::chrono::microseconds{dist(thread_data[id].rng)});
                }
            });
        });
}

#else

void work(thread_coordination::Context& ctx, PriorityQueue::Handle& handle) {
    ctx.execute_synchronized_blockwise_timed(
        "work", keys, keys + settings.num_threads * settings.num_operations,
        [&, id = ctx.get_id()](key_type* begin, key_type* end) {
            PriorityQueue::value_type retval;
            // Let retval escape
            asm volatile("" ::"g"(&retval));
            std::for_each(begin, end, [&, id](key_type k) {
                if (k == std::numeric_limits<std::uint32_t>::max()) {
                    while (!handle.try_extract_top(retval)) {
                        // Not expected due to prefill
                        ++thread_data[id].op_count.num_failed_deletions;
                    }
                    // "Use" memory to force write to retval
                    asm volatile("" ::: "memory");
                    ++thread_data[id].op_count.num_deletions;
                } else {
                    handle.push({k, k});
                    ++thread_data[id].op_count.num_insertions;
                }
            });
        });
}

#endif

struct Task {
    static void run(thread_coordination::Context ctx, PriorityQueue& pq) {
        if (ctx.is_main()) {
            std::clog << "Generating keys..." << std::flush;
        }
        generate_prefill_keys(ctx);
        generate_keys(ctx);

        if (ctx.is_main()) {
            std::clog << "done\nPrefilling..." << std::flush;
        }

        PriorityQueue::Handle handle = pq.get_handle();

        prefill(ctx, handle);

        if (ctx.is_main()) {
            std::clog << "done\nStarting the stress test..." << std::flush;
        }
        work(ctx, handle);
        if (ctx.is_main()) {
            std::clog << "done\n";
        }
    }

    static threading::thread_config get_config(
        thread_coordination::Context const& ctx) {
        threading::thread_config config;
        config.cpu_set.reset();
        config.cpu_set.set(ctx.get_id());
        return config;
    }
};

int main(int argc, char* argv[]) {
    std::clog << "Build configuration\n";
#ifndef NDEBUG
    std::clog << "DEBUG build\n";
#endif
#if defined THROUGHPUT_MODE
    std::clog << "Mode: Throughput\n";
#elif defined QUALITY_MODE
    std::clog << "Mode: Quality\n";
#ifdef USE_TSC
    std::clog << "Use TSC\n";
#endif
#endif
    std::clog << "L1 cache linesize (byte): " << L1_CACHE_LINESIZE << "\n\t";
    std::clog << "Pagesize (byte): " << PAGESIZE << "\n\t";
    std::clog << '\n';

    std::clog << "Command line: ";
    std::copy(argv, argv + argc,
              std::ostream_iterator<char const*>(std::clog, " "));
    std::clog << "\n\n";

    cxxopts::Options options(argv[0]);
    // clang-format off
    options.add_options()
      ("p,prefill", "Specify the number of elements to prefill the queue with "
       "(default: 1'000'000)", cxxopts::value<size_t>(), "NUMBER")
      ("n,ops", "Specify the number of operations per thread"
       "(default: 10'000'000)", cxxopts::value<std::size_t>(), "NUMBER")
      ("i,insert", "Specify the insert policy as one of \"uniform\", \"split\", \"producer\", \"alternating\" "
       "(default: uniform)", cxxopts::value<std::string>(), "ARG")
      ("j,threads", "Specify the number of threads "
       "(default: 4)", cxxopts::value<unsigned int>(), "NUMBER")
#ifdef QUALITY_MODE
      ("w,sleep", "Specify the sleep time between operations in microseconds"
       "(default: 0)", cxxopts::value<unsigned int>(), "NUMBER")
#endif
      ("d,distribution", "Specify the key distribution as one of \"uniform\", \"dijkstra\", \"ascending\", \"descending\", \"threadid\" "
       "(default: uniform)", cxxopts::value<std::string>(), "ARG")
      ("m,max", "Specify the max key "
       "(default: MAX)", cxxopts::value<key_type>(), "NUMBER")
      ("l,min", "Specify the min key "
       "(default: 0)", cxxopts::value<key_type>(), "NUMBER")
      ("s,seed", "Specify the initial seed"
       "(default: 0)", cxxopts::value<std::uint32_t>(), "NUMBER")
      ("c,factor", "The number of queues when using multiqueue or multififo"
       "(default: 4)", cxxopts::value<std::size_t>(), "NUMBER")
      ("k,stickiness", "The stickiness when using multiqueue or multififo supporting stickiness"
       "(default: 8)", cxxopts::value<unsigned int>(), "NUMBER")
      ("h,help", "Print this help");
    // clang-format on

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
        if (result.count("insert") > 0) {
            std::string policy = result["insert"].as<std::string>();
            if (policy == "uniform") {
                settings.insert_config.insert_policy = InsertPolicy::Uniform;
            } else if (policy == "alternating") {
                settings.insert_config.insert_policy =
                    InsertPolicy::Alternating;
            } else {
                std::cerr << "Unknown insert policy \"" << policy << "\"\n";
                return 1;
            }
        }
        if (result.count("distribution") > 0) {
            std::string dist = result["distribution"].as<std::string>();
            if (dist == "uniform") {
                settings.insert_config.key_distribution =
                    KeyDistribution::Uniform;
            } else if (dist == "ascending") {
                settings.insert_config.key_distribution =
                    KeyDistribution::Ascending;
            } else if (dist == "descending") {
                settings.insert_config.key_distribution =
                    KeyDistribution::Descending;
            } else if (dist == "dijkstra") {
                settings.insert_config.key_distribution =
                    KeyDistribution::Dijkstra;
            } else {
                std::cerr << "Unknown key distribution \"" << dist << "\"\n";
                return 1;
            }
        }
        if (result.count("threads") > 0) {
            settings.num_threads = result["threads"].as<unsigned int>();
        }
#ifdef QUALITY_MODE
        if (result.count("sleep") > 0) {
            settings.sleep_between_operations =
                std::chrono::microseconds{result["sleep"].as<unsigned int>()};
        }
#endif
        if (result.count("max") > 0) {
            settings.insert_config.max_key = result["max"].as<key_type>();
        }
        if (result.count("min") > 0) {
            settings.insert_config.min_key = result["min"].as<key_type>();
        }
        if (result.count("seed") > 0) {
            settings.seed = result["seed"].as<std::uint32_t>();
        }
        if (result.count("factor") > 0) {
            settings.pq_params.c = result["factor"].as<std::size_t>();
        }
        if (result.count("stickiness") > 0) {
            settings.pq_params.stickiness =
                result["stickiness"].as<unsigned int>();
        }
    } catch (cxxopts::OptionParseException const& e) {
        std::cerr << "Error parsing arguments: " << e.what() << '\n';
        std::cerr << options.help() << std::endl;
        return 1;
    }

    std::clog << "Settings: \n\t"
              << "Prefill: " << settings.prefill_size << "\n\t"
              << "Operations per thread: " << settings.num_operations << "\n\t"
              << "Threads: " << settings.num_threads << "\n\t"
              << "Insert policy: "
              << get_insert_policy_name(settings.insert_config.insert_policy)
              << "\n\t"
              << "Key distribution: "
              << get_key_distribution_name(
                     settings.insert_config.key_distribution)
              << "\n\t"
              << "Min key: " << settings.insert_config.min_key << "\n\t"
              << "Max key: " << settings.insert_config.max_key << "\n\t"
#ifdef QUALITY_MODE
              << "Sleep between operations: "
              << settings.sleep_between_operations.count() << " us\n\t"
#endif
              << "Dijkstra min increase: "
              << settings.insert_config.dijkstra_min_increase << "\n\t"
              << "Dijkstra max increase: "
              << settings.insert_config.dijkstra_max_increase << "\n\t"
              << "Seed: " << settings.seed;
    std::clog << "\n\n";

#ifdef QUALITY_MODE
    if (settings.num_threads > (1 << bits_for_thread_id) - 1) {
        std::cerr << "Too many threads, increase the number of thread bits!"
                  << std::endl;
        return 1;
    }
#endif
    xoroshiro256starstar rng;
    rng.seed(settings.seed);

    settings.pq_params.seed = rng();
    auto pq = util::create_pq<PriorityQueue>(
        settings.prefill_size, settings.num_threads, settings.pq_params);

    std::clog << "Using priority queue: " << pq.description() << "\n\n";

    thread_data = ::new ThreadData[settings.num_threads];
    for (std::size_t i = 0; i < settings.num_threads; ++i) {
        thread_data[i].rng.seed(rng());
        thread_data[i].inserter =
            InsertingStrategy<key_type>(settings.insert_config);
    }

    prefill_keys = ::new (std::align_val_t{L1_CACHE_LINESIZE})
        key_type[settings.prefill_size];

    keys = ::new (std::align_val_t{L1_CACHE_LINESIZE})
        key_type[settings.num_threads * settings.num_operations];
    thread_coordination::ThreadCoordinator coordinator{settings.num_threads};
    coordinator.run_task<Task>(std::ref(pq));
    coordinator.join();

    auto generate_prefill_duration =
        *coordinator.get_duration("generate_prefill_keys");
    auto generate_duration = *coordinator.get_duration("generate_keys");
    auto prefill_duration = *coordinator.get_duration("prefill");
    auto work_duration = *coordinator.get_duration("work");

    OperationCount sum_ops = std::transform_reduce(
        thread_data, thread_data + settings.num_threads, OperationCount{},
        std::plus<>(), [](auto const& d) { return d.op_count; });
    std::clog << "Insertions: " << sum_ops.num_insertions << '\n'
              << "Deletions: " << sum_ops.num_deletions << '\n'
              << "Failed deletions: " << sum_ops.num_failed_deletions << '\n';
    std::clog << "Generation time (s): " << std::setprecision(3)
              << std::chrono::duration<double>(generate_prefill_duration +
                                               generate_duration)
                     .count()
              << '\n';
    std::clog << "Prefill time (s): " << std::setprecision(3)
              << std::chrono::duration<double>(prefill_duration).count()
              << '\n';
    std::clog << "Work time (s): " << std::setprecision(3)
              << std::chrono::duration<double>(work_duration).count() << '\n';

#ifdef QUALITY_MODE
    std::cout << settings.num_threads << '\n';
    for (unsigned int t = 0; t < settings.num_threads; ++t) {
        for (auto const& [tick, key] : thread_data[t].ins_log) {
            std::cout << "i " << t << ' ' << tick << ' ' << key << '\n';
        }
    }

    for (unsigned int t = 0; t < settings.num_threads; ++t) {
        for (auto const& [tick, value] : thread_data[t].del_log) {
            std::cout << "d " << t << ' ' << tick << ' ' << get_thread_id(value)
                      << ' ' << get_elem_id(value) << '\n';
        }
    }

    for (unsigned int t = 0; t < settings.num_threads; ++t) {
        for (auto const& tick : thread_data[t].failed_del_log) {
            std::cout << "f " << t << ' ' << tick << '\n';
        }
    }

    std::cout << std::flush;

#else
    std::cout << "Ops/s: " << std::fixed << std::setprecision(2)
              << static_cast<double>(settings.num_threads *
                                     settings.num_operations) /
                     std::chrono::duration<double>(work_duration).count()
              << std::endl;
#endif
    ::operator delete[](keys, std::align_val_t{L1_CACHE_LINESIZE});
    ::operator delete[](prefill_keys, std::align_val_t{L1_CACHE_LINESIZE});
    delete[] thread_data;
    return 0;
}
