#include "utils/priority_queue_factory.hpp"
#include "utils/thread_coordination.hpp"
#include "utils/threading.hpp"

#include "external/cxxopts.hpp"

#include <time.h>
#include <x86intrin.h>
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
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

using key_type = unsigned long;
using value_type = unsigned long;

using PQFactory = util::PriorityQueueFactory<key_type, value_type, true>;
using PriorityQueue = typename PQFactory::type;

using Handle = typename PriorityQueue::handle_type;

static constexpr key_type PopOp =
    static_cast<value_type>(std::numeric_limits<std::uint32_t>::max());

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

struct Settings {
    std::size_t prefill_size = 1'000'000;
    std::size_t num_operations = 10'000'000;
    std::chrono::microseconds sleep_between_operations =
        std::chrono::microseconds::zero();
    unsigned int num_threads = 4;
    std::uint64_t seed = 1;
    util::PriorityQueueConfig pq_config;
    key_type min_key = 1;
    key_type max_key = PopOp >> 2;
    double pop_prob = 0.5;
};

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

struct PushLogEntry {
    tick_type tick;
    key_type key;
};

struct PopLogEntry {
    tick_type tick;
    value_type value;
};

struct alignas(2 * L1_CACHE_LINESIZE) ThreadData {
    unsigned int num_ops = 0;
    std::vector<PushLogEntry> push_log;
    std::vector<PopLogEntry> pop_log;
    std::vector<tick_type> failed_pop_log;
};

static Settings settings;
static ThreadData* thread_data;
static key_type** prefill_keys;
static key_type* keys;
static std::chrono::steady_clock::duration prefill_time;
static std::chrono::steady_clock::duration work_time;

bool check_operations() {
    auto pop_ops = static_cast<std::size_t>(
        std::count(keys, keys + settings.num_operations, PopOp));
    return settings.num_threads *
                   (settings.prefill_size / settings.num_threads) +
               settings.num_operations - pop_ops >=
           pop_ops;
}

struct Task {
    template <typename Generator>
    static void generate_prefill_keys(thread_coordination::Context& ctx,
                                      Generator& g) {
        auto dist = std::uniform_int_distribution<key_type>(settings.min_key,
                                                            settings.max_key);
        std::size_t prefill_by_thread =
            settings.prefill_size / settings.num_threads;
        prefill_keys[ctx.get_id()] = new key_type[prefill_by_thread];
        std::generate(prefill_keys[ctx.get_id()],
                      prefill_keys[ctx.get_id()] + prefill_by_thread,
                      [&]() { return dist(g); });
    }

    template <typename Generator>
    static void generate_keys(thread_coordination::Context& ctx, Generator& g) {
        auto key_dist = std::uniform_int_distribution<key_type>(
            settings.min_key, settings.max_key);
        auto pop_dist = std::uniform_real_distribution<>();
        ctx.execute_synchronized_blockwise(
            keys, keys + settings.num_operations,
            [&](key_type* begin, key_type* end) {
                std::generate(begin, end, [&]() {
                    return pop_dist(g) < settings.pop_prob ? PopOp
                                                           : key_dist(g);
                });
            });
    }

    static void prefill(thread_coordination::Context& ctx, Handle& handle) {
        ctx.execute_synchronized_timed(prefill_time, [id = ctx.get_id(),
                                                      &handle]() {
            std::size_t prefill_by_thread =
                settings.prefill_size / settings.num_threads;
            for (std::size_t i = 0; i < prefill_by_thread; ++i) {
                handle.push({prefill_keys[id][i], prefill_keys[id][i]});
                thread_data[id].push_log.push_back({0, prefill_keys[id][i]});
            }
        });
    }

    template <typename Generator>
    static void work(thread_coordination::Context& ctx, Handle& handle,
                     Generator& g) {
        unsigned int num_ops = 0;
        ctx.execute_synchronized_blockwise_timed(
            work_time, keys, keys + settings.num_operations,
            [&handle, id = ctx.get_id(), &num_ops, &g](key_type* begin,
                                                       key_type* end) {
                for (auto it = begin; it != end; ++it) {
                    if (*it == PopOp) {
                        PriorityQueue::value_type retval;
                        while (!handle.try_pop(retval)) {
                            auto tick = get_tick();
                            thread_data[id].failed_pop_log.push_back(tick);
                        }
                        auto tick = get_tick();
                        thread_data[id].pop_log.push_back(
                            {tick, retval.second});
                    } else {
                        value_type value =
                            to_value(id, thread_data[id].push_log.size());
                        handle.push({*it, value});
                        auto tick = get_tick();
                        thread_data[id].push_log.push_back({tick, *it});
                    }
                    if (settings.sleep_between_operations !=
                        std::chrono::microseconds::zero()) {
                        auto dist =
                            std::uniform_int_distribution<std::uint64_t>(
                                1,
                                static_cast<std::uint64_t>(
                                    settings.sleep_between_operations.count()));
                        std::this_thread::sleep_for(
                            std::chrono::microseconds{dist(g)});
                    }
                }
                num_ops += static_cast<unsigned int>(end - begin);
            });
        thread_data[ctx.get_id()].num_ops = num_ops;
    }

    static void run(thread_coordination::Context ctx, std::uint64_t* seeds,
                    PriorityQueue& pq) {
        if (ctx.get_id() == 0) {
            std::clog << "Generating keys..." << std::flush;
        }

        std::default_random_engine rng(seeds[ctx.get_id()]);
        generate_prefill_keys(ctx, rng);
        generate_keys(ctx, rng);

        if (ctx.get_id() == 0) {
            std::clog << "done" << std::endl;
        }

        Handle handle =
            ctx.execute_exclusive([&pq]() { return pq.get_handle(); });

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

        work(ctx, handle, rng);

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
    std::filesystem::path log_file;

    cxxopts::Options options(argv[0]);
    options.add_options()
        // clang-format off
        ("h,help", "Print this help")
        ("p,prefill", "Specify the number of elements to prefill the queue with "
         "(default: 1'000'000)", cxxopts::value<size_t>(), "NUMBER")
        ("n,ops", "Specify the number of operations"
         "(default: 10'000'000)", cxxopts::value<std::size_t>(), "NUMBER")
        ("j,threads", "Specify the number of threads "
         "(default: 4)", cxxopts::value<unsigned int>(), "NUMBER")
        ("d,pop-prob", "Specify the probability of pops"
         "(default: 0.5)", cxxopts::value<double>(), "NUMBER")
        ("w,sleep", "Specify the sleep time between operations in microseconds "
         "(default: 0)", cxxopts::value<unsigned int>(), "NUMBER")
        ("m,max", "Specify the max key "
         "(default: MAX)", cxxopts::value<key_type>(), "NUMBER")
        ("l,min", "Specify the min key "
         "(default: 0)", cxxopts::value<key_type>(), "NUMBER")
        ("s,seed", "Specify the initial seed "
         "(default: 0)", cxxopts::value<std::uint32_t>(), "NUMBER")
        ("c,factor", "The number of queues when using multiqueue or multififo"
         "(default: 4)", cxxopts::value<std::size_t>(), "NUMBER")
        ("k,stickiness", "The stickiness (only has effect for some multiqueue variants)"
         "(default: 8)", cxxopts::value<unsigned int>(), "NUMBER")
        ("o,output", "Output data in csv (comma-separated)", cxxopts::value<std::filesystem::path>(out_file), "PATH")
        ("f,log-file", "Path to the logfile", cxxopts::value<std::filesystem::path>(log_file)->default_value("log.txt"), "PATH")
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
        if (result.count("sleep") > 0) {
            settings.sleep_between_operations =
                std::chrono::microseconds{result["sleep"].as<unsigned int>()};
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
            settings.pq_config.stickiness =
                result["stickiness"].as<unsigned int>();
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
              << "Sleep between operations: "
              << settings.sleep_between_operations.count() << " us\n\t"
              << "Pop probability: " << std::fixed << settings.pop_prob << '\n'
              << "Min key: " << settings.min_key << '\n'
              << "Max key: " << settings.max_key << '\n'
              << "Seed: " << settings.seed << '\n';
    std::cout << '\n';

    if (settings.num_threads > (1 << bits_for_thread_id) - 1) {
        std::cerr << "Too many threads!" << std::endl;
        return 1;
    }
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
    thread_coordination::run_task<Task>(settings.num_threads,
                                        thread_seeds.get(), std::ref(pq));
    unsigned int total_failed_pops =
        std::accumulate(thread_data, thread_data + settings.num_threads, 0u,
                        [](unsigned int sum, auto const& d) {
                            return sum + d.failed_pop_log.size();
                        });
    auto [min_thread, max_thread] =
        std::minmax_element(thread_data, thread_data + settings.num_threads,
                            [](auto const& lhs, auto const& rhs) {
                                return lhs.num_ops < rhs.num_ops;
                            });
    std::cout << "Least/Most operations by any thread: " << min_thread->num_ops
              << '/' << max_thread->num_ops << '\n';
    std::cout << "Failed pops: " << total_failed_pops << '\n';
    std::cout << "Prefill time (s): " << std::setprecision(3)
              << std::chrono::duration<double>(prefill_time).count() << '\n';
    std::cout << "Work time (s): " << std::setprecision(3)
              << std::chrono::duration<double>(work_time).count() << '\n';

    std::clog << "Writing log file..." << std::flush;
    auto log_out = std::ofstream{log_file};
    if (!log_out.is_open()) {
        std::cerr << "Could not open log file\n";
        return 1;
    }
    log_out << settings.num_threads << '\n';
    for (unsigned int t = 0; t < settings.num_threads; ++t) {
        for (auto const& [tick, key] : thread_data[t].push_log) {
            log_out << "i " << t << ' ' << tick << ' ' << key << '\n';
        }
    }

    for (unsigned int t = 0; t < settings.num_threads; ++t) {
        for (auto const& [tick, value] : thread_data[t].pop_log) {
            log_out << "d " << t << ' ' << tick << ' ' << get_thread_id(value)
                    << ' ' << get_elem_id(value) << '\n';
        }
    }

    for (unsigned int t = 0; t < settings.num_threads; ++t) {
        for (auto const& tick : thread_data[t].failed_pop_log) {
            log_out << "f " << t << ' ' << tick << '\n';
        }
    }

    log_out << std::flush;
    log_out.close();
    std::clog << "done" << std::endl;

    if (!out_file.empty()) {
        auto out = std::ofstream{out_file};
        if (!out.is_open()) {
            std::cerr << "Could not open file to write out results\n";
            return 1;
        }
        out << "prefill,ops,threads,pop_prob,min_key,max_key,"
               "seed,min_thread_ops,max_thread_ops,failed_pops,prefill_time,"
               "work_time"
            << '\n';
        out << settings.prefill_size << ',' << settings.num_operations << ','
            << settings.num_threads << ',' << settings.pop_prob << ','
            << settings.min_key << ',' << settings.max_key << ','
            << settings.seed << ',' << min_thread->num_ops << ','
            << max_thread->num_ops << ',' << total_failed_pops << ','
            << std::fixed << std::setprecision(2)
            << std::chrono::duration<double>(prefill_time).count() << ','
            << std::fixed << std::setprecision(2)
            << std::chrono::duration<double>(work_time).count()
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
