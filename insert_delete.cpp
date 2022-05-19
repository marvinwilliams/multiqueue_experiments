#include "utils/priority_queue_factory.hpp"
#include "utils/thread_coordination.hpp"
#include "utils/threading.hpp"
#include "utils/xoroshiro256starstar.hpp"

#include "external/cxxopts.hpp"

#include <time.h>
#ifdef __x86_64__
#include <x86intrin.h>
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
#include <numeric>
#include <optional>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

#ifndef L1_CACHE_LINESIZE
#define L1_CACHE_LINESIZE 64
#endif

using key_type = unsigned long;
using value_type = unsigned long;

using PriorityQueue =
    typename util::PriorityQueueFactory<key_type, value_type>::type;

struct Settings {
    std::size_t prefill_size = 1'000'000;
    std::size_t num_operations = 10'000'000;
    unsigned int num_threads = 4;
    std::uint64_t seed = 1;
    util::PriorityQueueParameters pq_params;
    key_type min_key = 1;
    key_type max_key = (1ul << 32) - 2;
    bool json_output = false;
};

struct alignas(L1_CACHE_LINESIZE) ThreadData {
    xoroshiro256starstar rng;
};

static Settings settings;
alignas(L1_CACHE_LINESIZE) static ThreadData* thread_data;
alignas(L1_CACHE_LINESIZE) static key_type* prefill_keys;
alignas(L1_CACHE_LINESIZE) static key_type* keys;

std::atomic<std::size_t> counter = 2;

void generate_keys(thread_coordination::Context& ctx) {
    ctx.execute_synchronized_blockwise_timed(
        "generate_prefill_keys", prefill_keys,
        prefill_keys + settings.prefill_size,
        [id = ctx.get_id()](key_type* begin, key_type* end) {
            std::generate(begin, end, [id]() {
                /* return 2; */
                /* return counter++; */
                auto dist = std::uniform_int_distribution<key_type>(
                    settings.min_key, settings.max_key);
                return dist(thread_data[id].rng);
            });
        });
    ctx.execute_synchronized_blockwise_timed(
        "generate_keys", keys,
        keys + settings.num_threads * settings.num_operations,
        [id = ctx.get_id()](key_type* begin, key_type* end) {
            std::generate(begin, end, [id]() {
                /* return 2; */
                /* return counter++; */
                auto dist = std::uniform_int_distribution<key_type>(
                    settings.min_key, settings.max_key);
                return dist(thread_data[id].rng);
            });
        });
}

void prefill(thread_coordination::Context& ctx, PriorityQueue::Handle& handle) {
    ctx.execute_synchronized_blockwise_timed(
        "prefill", prefill_keys, prefill_keys + settings.prefill_size,
        [&handle](key_type* begin, key_type* end) {
            std::for_each(begin, end, [&handle](key_type k) {
                handle.push({k, k});
            });
        });
}

void insert_all(thread_coordination::Context& ctx,
                PriorityQueue::Handle& handle) {
    ctx.execute_synchronized_blockwise_timed(
        "insert", keys, keys + settings.num_threads * settings.num_operations,
        [&handle, &ctx](key_type* begin, key_type* end) {
            /* ctx.write(std::cout) << end - begin << '\n'; */
            std::for_each(begin, end, [&handle](key_type k) {
                handle.push({k, k});
            });
        });
}

void delete_all(thread_coordination::Context& ctx,
                PriorityQueue::Handle& handle) {
    static constexpr std::size_t block_size = 4096;
    static std::atomic_size_t index{0};

    ctx.execute_synchronized_timed(
        "delete",
        [&handle, &ctx](std::size_t n) {
            while (index.load(std::memory_order_relaxed) < n) {
                std::size_t i = 0;
                for (; i < block_size; ++i) {
                    PriorityQueue::value_type retval;
                    // Let retval escape
                    asm volatile("" ::"g"(&retval));
                    if (!handle.try_extract_top(retval)) {
                        break;
                    }
                    // "Use" memory to force write to retval
                    asm volatile("" ::: "memory");
                }
                /* ctx.write(std::cout) << "s " << num_delete << '\n'; */
                index.fetch_add(i, std::memory_order_relaxed);
            }
        },
        settings.num_threads * settings.num_operations);
}

struct Task {
    static void run(thread_coordination::Context ctx, PriorityQueue& pq) {
        if (ctx.is_main()) {
            std::clog << "Generating keys..." << std::endl;
        }

        generate_keys(ctx);

        PriorityQueue::Handle handle = pq.get_handle();

        if (ctx.is_main()) {
            std::clog << "Prefilling..." << std::endl;
        }

        prefill(ctx, handle);

        if (ctx.is_main()) {
            std::clog << "Inserting..." << std::endl;
        }

        insert_all(ctx, handle);

        if (ctx.is_main()) {
            std::clog << "Deleting..." << std::endl;
        }

        delete_all(ctx, handle);

        if (ctx.is_main()) {
            std::clog << std::endl;
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
#ifndef NDEBUG
    std::clog << "DEBUG build\n";
#endif
    std::clog << "Command line: ";
    std::copy(argv, argv + argc,
              std::ostream_iterator<char const*>(std::clog, " "));
    std::clog << '\n' << '\n';

    cxxopts::Options options(argv[0]);
    // clang-format off
    options.add_options()
      ("p,prefill", "Specify the number of elements to prefill the queue with "
       "(default: 1'000'000)", cxxopts::value<size_t>(), "NUMBER")
      ("n,ops", "Specify the number of operations per thread"
       "(default: 10'000'000)", cxxopts::value<std::size_t>(), "NUMBER")
      ("j,threads", "Specify the number of threads "
       "(default: 4)", cxxopts::value<unsigned int>(), "NUMBER")
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
      ("o,json", "Produce json output")
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
        if (result.count("threads") > 0) {
            settings.num_threads = result["threads"].as<unsigned int>();
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
            settings.pq_params.c = result["factor"].as<std::size_t>();
        }
        if (result.count("stickiness") > 0) {
            settings.pq_params.stickiness =
                result["stickiness"].as<unsigned int>();
        }
        if (result.count("json") > 0) {
            settings.json_output = true;
        }
    } catch (cxxopts::OptionParseException const& e) {
        std::cerr << "Error parsing arguments: " << e.what() << '\n';
        std::cerr << options.help() << std::endl;
        return 1;
    }

    if (settings.json_output) {
        std::cout << "{\"settings\": {"
                  << "\"prefill\": " << settings.prefill_size << ','
                  << "\"ops_per_thread\": " << settings.num_operations << ','
                  << "\"threads\": " << settings.num_threads << ','
                  << "\"min_key\": " << settings.min_key << ','
                  << "\"max_key\": " << settings.max_key << ','
                  << "\"seed\": " << settings.seed;
        std::cout << "}";
    } else {
        std::cout << "Settings\n";
        std::cout << "  Prefill: " << settings.prefill_size << '\n'
                  << "  Ops per thread: " << settings.num_operations << '\n'
                  << "  Threads: " << settings.num_threads << '\n'
                  << "  Min key: " << settings.min_key << '\n'
                  << "  Max key: " << settings.max_key << '\n'
                  << "  Seed: " << settings.seed << '\n';
        std::cout << '\n';
    }
    xoroshiro256starstar rng;
    rng.seed(settings.seed);

    settings.pq_params.seed = rng();
    auto pq = util::create_pq<PriorityQueue>(
        settings.prefill_size, settings.num_threads, settings.pq_params);

    if (settings.json_output) {
    } else {
        std::cout << "Priority queue description\n  " << pq.description()
                  << '\n'
                  << '\n';
    }

    thread_data = ::new ThreadData[settings.num_threads];
    for (std::size_t i = 0; i < settings.num_threads; ++i) {
        thread_data[i].rng.seed(rng());
    }

    prefill_keys = ::new (std::align_val_t{L1_CACHE_LINESIZE})
        key_type[settings.prefill_size];

    keys = ::new (std::align_val_t{L1_CACHE_LINESIZE})
        key_type[settings.num_threads * settings.num_operations];
    thread_coordination::ThreadCoordinator coordinator{settings.num_threads};
    coordinator.run_task<Task>(std::ref(pq));
    coordinator.join();

    ::delete[] keys;
    ::delete[] prefill_keys;
    ::delete[] thread_data;

    auto generate_prefill_duration =
        *coordinator.get_duration("generate_prefill_keys");
    auto generate_duration = *coordinator.get_duration("generate_keys");
    auto prefill_duration = *coordinator.get_duration("prefill");
    auto insert_duration = *coordinator.get_duration("insert");
    auto delete_duration = *coordinator.get_duration("delete");

    if (settings.json_output) {
        std::cout << ", \"results\": {";
        std::cout
            << std::fixed << "\"gen_prefill_keys_time\": "
            << std::chrono::duration<double>(generate_prefill_duration).count()
            << ",\"gen_workload_keys_time\": "
            << std::chrono::duration<double>(generate_duration).count()
            << ",\"prefill_time\": "
            << std::chrono::duration<double>(prefill_duration).count()
            << ",\"insert_time\": "
            << std::chrono::duration<double>(insert_duration).count()
            << ",\"delete_time\": "
            << std::chrono::duration<double>(delete_duration).count()
            << ",\"insert_throughput\": "
            << static_cast<double>(settings.num_operations) /
                   std::chrono::duration<double>(insert_duration).count()
            << ",\"delete_throughput\": "
            << static_cast<double>(settings.num_operations) /
                   std::chrono::duration<double>(delete_duration).count();
        std::cout << "}}\n";
    } else {
        std::cout
            << std::fixed << std::setprecision(3)
            << "Prefill key generation time (s): "
            << std::chrono::duration<double>(generate_prefill_duration).count()
            << '\n'
            << "Workload key generation time (s): "
            << std::chrono::duration<double>(generate_duration).count() << '\n'
            << "Prefill time (s): "
            << std::chrono::duration<double>(prefill_duration).count() << '\n'
            << "Insert time (s): "
            << std::chrono::duration<double>(insert_duration).count() << '\n'

            << "Delete time (s): "
            << std::chrono::duration<double>(delete_duration).count() << '\n'

            << "Insert throughput (ops/t/s): "
            << static_cast<double>(settings.num_operations) /
                   std::chrono::duration<double>(insert_duration).count()
            << '\n'

            << "Delete throughput (ops/t/s): "
            << static_cast<double>(settings.num_operations) /
                   std::chrono::duration<double>(delete_duration).count()
            << '\n';
    }
    return 0;
}
