#include "priority_queue_factory.hpp"
#include "thread_coordination.hpp"

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

#ifdef WITH_PAPI
extern "C" {
#include <papi.h>
#include <pthread.h>
}
#endif

using key_type = DefaultMinPriorityQueue::key_type;
using Handle = DefaultMinPriorityQueue::Handle;

#ifdef WITH_PAPI
// These are the event names used for the papi library to measure cache misses
// As these are hardware specific, you need to modify them or use generic PAPI events
static auto const L1d_cache_miss_event_name = "perf_raw::rc860";
static auto const L2_cache_miss_event_name = "perf_raw::r0864";

bool start_performance_counter(int& event_set) {
    if (int ret = PAPI_register_thread(); ret != PAPI_OK) {
        return false;
    }
    if (int ret = PAPI_create_eventset(&event_set); ret != PAPI_OK) {
        return false;
    }
    int code{};
    if (int ret = PAPI_event_name_to_code(L1d_cache_miss_event_name, &code); ret != PAPI_OK) {
        return false;
    }
    if (int ret = PAPI_add_event(event_set, code); ret != PAPI_OK) {
        return false;
    }
    if (int ret = PAPI_event_name_to_code(L2_cache_miss_event_name, &code); ret != PAPI_OK) {
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

struct Settings {
    enum class WorkMode { Mixed, Split };
    enum class ElementDistribution { Uniform, Ascending, Descending };

    auto work_mode_str() const {
        switch (work_mode) {
            case WorkMode::Mixed:
                return "mixed";
            case WorkMode::Split:
                return "split";
        }
        return "unknown";
    }
    auto element_distribution_str() const {
        switch (element_distribution) {
            case ElementDistribution::Uniform:
                return "uniform";
            case ElementDistribution::Ascending:
                return "ascending";
            case ElementDistribution::Descending:
                return "descending";
        }
        return "unknown";
    }

    int num_threads = 4;
    std::size_t prefill_per_thread = 1 << 20;
    std::size_t operations_per_thread = 1 << 24;
    key_type min_key = 1;
    key_type max_key = key_type{1} << 30;
    double pop_probability = 0.5;
    int num_push_threads = 1;
    WorkMode work_mode = WorkMode::Mixed;
    ElementDistribution element_distribution = ElementDistribution::Uniform;
#ifdef WITH_PAPI
    bool enable_performance_counter = false;
#endif
    int seed = 1;
};

struct Result {
    std::atomic<thread_coordination::duration_type> work_time{};
    std::atomic<long long> num_failed_pops{0};
    std::atomic<long long> num_pops{0};
#ifdef WITH_PAPI
    std::atomic<long long> l1d_cache_misses{0};
    std::atomic<long long> l2_cache_misses{0};
#endif
#ifdef MQ_COUNT_STATS
    std::atomic<long long> num_locking_failed{0};
    std::atomic<long long> num_resets{0};
    std::atomic<long long> use_counts{0};
#endif

    void update_work_time(thread_coordination::duration_type t) {
        auto old = work_time.load(std::memory_order_relaxed);
        while (t > old && !work_time.compare_exchange_weak(old, t, std::memory_order_relaxed)) {
        }
    }
};

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
    [[nodiscard]] constexpr key_type get_key() const noexcept {
        assert(!is_pop());
        return data;
    }
};

template <typename OpIt>
void execute_mixed(thread_coordination::Context const& ctx, Handle& handle, OpIt begin, OpIt end, Result& result) {
    long long num_failed_pops = 0;
    volatile DefaultMinPriorityQueue::mapped_type val;
    DefaultMinPriorityQueue::value_type retval;
    auto t = ctx.execute_synchronized_blockwise(
        begin, end, [&handle, &num_failed_pops, &val, &retval](auto block_begin, auto block_end) {
            for (auto it = block_begin; it != block_end; ++it) {
                if (it->is_pop()) {
                    while (!handle.try_pop(retval)) {
                        ++num_failed_pops;
                    }
                    val = retval.second;
                } else {
                    key_type key = it->get_key();
                    handle.push({key, key});
                }
            }
        });
    result.num_failed_pops += num_failed_pops;
    result.update_work_time(t);
}

template <typename OpIt>
void execute_split_push(thread_coordination::Context const& ctx, Handle& handle, OpIt begin, OpIt end, Result& result) {
    long long num_failed_pops = 0;
    auto t =
        ctx.execute_synchronized_blockwise(begin, end, [&handle, &num_failed_pops](auto block_begin, auto block_end) {
            for (auto it = block_begin; it != block_end; ++it) {
                key_type key = it->get_key();
                handle.push({key, key});
            }
        });
    result.update_work_time(t);
}

void execute_split_pop(thread_coordination::Context const& ctx, Handle& handle, Result& result,
                       std::size_t num_elements) {
    long long num_failed_pops = 0;
    auto t = ctx.execute_synchronized([&handle, &result, &num_failed_pops, num_elements]() {
        volatile DefaultMinPriorityQueue::mapped_type val;
        DefaultMinPriorityQueue::value_type retval;
        long long num_pops;
        do {
            num_pops = 0;
            while (handle.try_pop(retval)) {
                val = retval.second;
                ++num_pops;
            }
            ++num_failed_pops;
        } while (result.num_pops.fetch_add(num_pops, std::memory_order_relaxed) + num_pops <
                 static_cast<long long>(num_elements));
    });
    assert(result.num_pops.load(std::memory_order_relaxed) == static_cast<long long>(num_elements));
    result.update_work_time(t);
    result.num_failed_pops += num_failed_pops;
}

void benchmark_thread(thread_coordination::Context ctx, Settings const& settings, DefaultMinPriorityQueue& pq,
                      std::vector<Operation> operations, Result& result) {
    ctx.synchronize([]() { std::clog << "Generating operations..." << std::flush; });

    std::seed_seq seed{settings.seed, ctx.get_id()};
    std::default_random_engine rng(seed);

    auto key_dist = std::uniform_int_distribution<key_type>(settings.min_key, settings.max_key);
    auto pop_dist = std::bernoulli_distribution(settings.pop_probability);

    // Generate workload
    auto start_index = static_cast<std::size_t>(ctx.get_id()) * settings.operations_per_thread;
    if (settings.element_distribution == Settings::ElementDistribution::Uniform) {
        std::generate_n(operations.begin() + static_cast<std::ptrdiff_t>(start_index), settings.operations_per_thread,
                        [&key_dist, &pop_dist, &rng, &settings]() {
                            return settings.work_mode == Settings::WorkMode::Mixed && pop_dist(rng)
                                ? Operation{}
                                : Operation{key_dist(rng)};
                        });
    } else {
        for (auto i = start_index; i < start_index + settings.operations_per_thread; ++i) {
            operations[i] = settings.work_mode == Settings::WorkMode::Mixed && pop_dist(rng)
                ? Operation{}
                : Operation{settings.min_key +
                            ((settings.element_distribution == Settings::ElementDistribution::Ascending
                                  ? i
                                  : operations.size() - i - 1) *
                             (settings.max_key - settings.min_key)) /
                                operations.size()};
        }
    }

    ctx.synchronize([]() { std::clog << "done\n"; });

    Handle handle = pq.get_handle(ctx.get_id());

    if (settings.work_mode == Settings::WorkMode::Mixed && settings.prefill_per_thread > 0) {
        if (ctx.get_id() == 0) {
            std::clog << "Prefilling..." << std::flush;
        }
        ctx.execute_synchronized([&handle, &key_dist, &rng, &settings]() {
            for (std::size_t i = 0; i < settings.prefill_per_thread; ++i) {
                auto key = key_dist(rng);
                handle.push({key, key});
            }
        });
        if (ctx.get_id() == 0) {
            std::clog << "done\n";
        }
    }
    if (ctx.get_id() == 0) {
        std::clog << "Working..." << std::flush;
    }

#ifdef MQ_COUNT_STATS
    handle.stats.num_locking_failed = 0;
    handle.stats.num_resets = 0;
    handle.stats.use_counts = 0;
#endif
#ifdef WITH_PAPI
    bool papi_started = false;
    int event_set = PAPI_NULL;
    if (settings.enable_performance_counter) {
        papi_started = start_performance_counter(event_set);
        if (!papi_started) {
            ctx.write(std::cerr) << "Failed to start counters\n";
        }
    }
#endif
    if (settings.work_mode == Settings::WorkMode::Mixed) {
        execute_mixed(ctx, handle, operations.begin(), operations.end(), result);
    } else {
        if (ctx.get_id() < settings.num_push_threads) {
            execute_split_push(ctx, handle, operations.begin(), operations.end(), result);
        } else {
            execute_split_pop(ctx, handle, result, operations.size());
        }
    }
#ifdef WITH_PAPI
    if (papi_started) {
        long long counter[2]{};
        if (int ret = PAPI_stop(event_set, counter); ret != PAPI_OK) {
            ctx.write(std::cerr) << "Failed to stop counters\n";
        } else {
            data.l1d_cache_misses += counter[0];
            data.l2_cache_misses += counter[1];
        }
    }
#endif
#ifdef MQ_COUNT_STATS
    data.num_locking_failed += handle.stats.num_locking_failed;
    data.num_resets += handle.stats.num_resets;
    data.use_counts += handle.stats.use_counts;
#endif
    if (ctx.get_id() == 0) {
        std::clog << "done\n" << std::endl;
    }
}

int main(int argc, char* argv[]) {
#ifndef NDEBUG
    std::clog << "Build type: Debug\n";
#else
    std::clog << "Build type: Release\n";
#endif
#ifdef WITH_PAPI
    std::clog << "Performance counter: enabled\n";
#else
    std::clog << "Performance counter: disabled\n";
#endif
    std::clog << "L1 cache linesize (bytes): " << L1_CACHE_LINESIZE << '\n';
    std::clog << "Pagesize (bytes): " << PAGESIZE << '\n';
    std::clog << '\n';

    std::clog << "Command line:";
    for (int i = 0; i < argc; ++i) {
        std::clog << ' ' << argv[i];
    }
    std::clog << '\n';

    Settings settings{};
    cxxopts::Options options(argv[0]);
    // clang-format off
    options.add_options()
        ("h,help", "Print this help")
        ("j,threads", "The number of threads", cxxopts::value<int>(settings.num_threads), "NUMBER")
        ("p,prefill", "The prefill per thread", cxxopts::value<std::size_t>(settings.prefill_per_thread), "NUMBER")
        ("n,ops", "The number of operations per thread", cxxopts::value<std::size_t>(settings.operations_per_thread), "NUMBER")
        ("l,min", "Specify the min key", cxxopts::value<key_type>(settings.min_key), "NUMBER")
        ("m,max", "Specify the max key", cxxopts::value<key_type>(settings.max_key), "NUMBER")
        ("o,work-mode", "Specify the work mode ([m]ixed, [s]plit)", cxxopts::value<char>(), "STRING")
        ("e,element-distribution", "Specify the element distribution ([u]niform, [a]scending, [d]escending)", cxxopts::value<char>(), "STRING")
        ("d,pop-prob", "Specify the probability of pops in [0,1] in mixed mode", cxxopts::value<double>(settings.pop_probability), "NUMBER")
        ("i,push-threads", "The number of pushing threads in split mode", cxxopts::value<int>(settings.num_push_threads), "NUMBER")
        ("s,seed", "Specify the initial seed", cxxopts::value<int>(settings.seed), "NUMBER")
#ifdef WITH_PAPI
        ("r,pc", "Enable Performance counter", cxxopts::value<>(settings.enable_performance_counter))
#endif
        ;
    // clang-format on
    add_pq_options(options);
    cxxopts::ParseResult args;
    try {
        args = options.parse(argc, argv);
    } catch (cxxopts::OptionParseException const& e) {
        std::cerr << "Error parsing arguments: " << e.what() << '\n';
        std::cerr << options.help() << std::endl;
        return 1;
    }
    if (args.count("help") > 0) {
        std::cerr << options.help() << std::endl;
        return 0;
    }
    if (args.count("work-mode") > 0) {
        if (args["work-mode"].as<char>() == 'm') {
            settings.work_mode = Settings::WorkMode::Mixed;
        } else if (args["work-mode"].as<char>() == 's') {
            settings.work_mode = Settings::WorkMode::Split;
        } else {
            std::cerr << "Invalid work mode: " << args["work-mode"].as<char>() << '\n';
            std::cerr << options.help() << std::endl;
            return 1;
        }
    }
    if (args.count("element-distribution") > 0) {
        if (args["element-distribution"].as<char>() == 'u') {
            settings.element_distribution = Settings::ElementDistribution::Uniform;
        } else if (args["element-distribution"].as<char>() == 'a') {
            settings.element_distribution = Settings::ElementDistribution::Ascending;
        } else if (args["element-distribution"].as<char>() == 'd') {
            settings.element_distribution = Settings::ElementDistribution::Descending;
        } else {
            std::cerr << "Invalid element distribution: " << args["element-distribution"].as<char>() << '\n';
            std::cerr << options.help() << std::endl;
            return 1;
        }
    }

    std::clog << "Threads: " << settings.num_threads << '\n'
              << "Prefill per thread: " << settings.prefill_per_thread << '\n'
              << "Operations per thread: " << settings.operations_per_thread << '\n'
              << "Min key: " << settings.min_key << '\n'
              << "Max key: " << settings.max_key << '\n'
              << "Pop probability: " << std::fixed << std::setprecision(2) << settings.pop_probability << '\n'
              << "Push threads: " << settings.num_push_threads << '\n'
              << "Operation mode: " << settings.work_mode_str() << '\n'
              << "Element distribution: " << settings.element_distribution_str() << '\n'
              << "Seed: " << settings.seed << '\n'
              << '\n';

#ifdef WITH_PAPI
    if (settings.enable_performance_counter) {
        if (int ret = PAPI_library_init(PAPI_VER_CURRENT); ret != PAPI_VER_CURRENT) {
            std::cerr << "Error initializing PAPI\n";
            return 1;
        }
        if (int ret = PAPI_thread_init((unsigned long (*)())(pthread_self)); ret != PAPI_OK) {
            std::cerr << "Error initializing threads for PAPI\n";
            return 1;
        }
        if (int ret = PAPI_query_named_event(L1d_cache_miss_event_name); ret != PAPI_OK) {
            std::cerr << "Cannot measure event '" << L1d_cache_miss_event_name << "'\n";
            return 1;
        }
        if (int ret = PAPI_query_named_event(L2_cache_miss_event_name); ret != PAPI_OK) {
            std::cerr << "Cannot measure event '" << L2_cache_miss_event_name << "'\n";
            return 1;
        }
    }
#endif

// TODO: Fix capacity
    auto pq = create_pq<DefaultMinPriorityQueue>(
        settings.num_threads, settings.prefill_per_thread * static_cast<std::size_t>(settings.num_threads), args);

    std::clog << "Priority queue: " << pq_name << '\n';

    Result result;
    std::vector<Operation> operations(static_cast<std::size_t>(settings.num_threads) * settings.operations_per_thread);
    thread_coordination::TaskHandle task_handle(settings.num_threads, benchmark_thread, settings, std::ref(pq),
                                                std::ref(operations), std::ref(result));
    task_handle.wait();

    std::clog << "Work time (s): " << std::setprecision(3)
              << std::chrono::duration<double>(result.work_time.load()).count() << '\n';
    std::clog << "Failed pops: " << result.num_failed_pops << '\n';
#ifdef WITH_PAPI
    if (settings.enable_performance_counter) {
        std::clog << "L1d cache misses: " << result.l1d_cache_misses << '\n';
        std::clog << "L2 cache misses: " << result.l2_cache_misses << '\n';
    }
#endif
#ifdef MQ_COUNT_STATS
    std::clog << "Failed locks per operation: "
              << static_cast<double>(result.num_locking_failed) /
            static_cast<double>(settings.num_threads * settings.operations_per_thread)
              << '\n';
    std::clog << "Average queue use count: "
              << static_cast<double>(result.use_counts) / static_cast<double>(result.num_resets) << '\n';
#endif

    std::cout << settings.num_threads << ',' << settings.prefill_per_thread << ',' << settings.operations_per_thread
              << ',' << settings.min_key << ',' << settings.max_key << ',' << settings.pop_probability << ','
              << settings.num_push_threads << ',' << settings.work_mode_str() << ','
              << settings.element_distribution_str() << ',' << settings.seed << ',' << std::fixed
              << std::setprecision(3) << std::chrono::duration<double>(result.work_time.load()).count() << ','
              << result.num_failed_pops;
#ifdef WITH_PAPI
    if (settings.enable_performance_counter) {
        std::cout << ',' << benchmark_data.l1d_cache_misses << ',' << benchmark_data.l2_cache_misses;
    } else {
        std::cout << ",n/a,n/a";
    }
#else
    std::cout << ",n/a,n/a";
#endif
#ifdef MQ_COUNT_STATS
    std::cout << ',' << benchmark_data.num_resets << ',' << benchmark_data.use_counts;
#else
    std::cout << ",n/a,n/a";
#endif
    std::cout << std::endl;
    return 0;
}
