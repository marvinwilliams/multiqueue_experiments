#include "priority_queue_factory.hpp"
#include "thread_coordination.hpp"
#ifdef QUALITY
#include "evaluate_quality.hpp"
#endif

#include "cxxopts.hpp"

#ifdef WITH_PAPI
extern "C" {
#include <papi.h>
#include <pthread.h>
}
#endif

#include <atomic>
#include <cassert>
#include <chrono>
#include <cstddef>
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
#ifdef QUALITY
#include <filesystem>
#endif

#ifdef USE_UINT8
#ifdef QUALITY
#error "QUALITY does not work with UINT8"
#endif
#ifndef PQ_HAS_GENERIC
#error "UINT8 is defined, but PQ is not generic"
#endif
using PriorityQueue = GenericMinPriorityQueue<std::uint8_t, std::uint8_t>;
#else
using PriorityQueue = DefaultMinPriorityQueue;
#endif
using key_type = PriorityQueue::key_type;
using handle_type = PriorityQueue::handle_type;

using thread_coordination::timepoint_type;

#ifdef QUALITY
std::uint64_t get_tick() noexcept {
    timespec ts;
    clock_gettime(CLOCK_MONOTONIC_RAW, &ts);
    return static_cast<std::uint64_t>(ts.tv_sec) * 1000000000 + static_cast<std::uint64_t>(ts.tv_nsec);
}
#endif

#ifdef WITH_PAPI
// These are the event names used for the papi library to measure cache misses
// As these are hardware specific, you need to modify them or use generic PAPI events
static auto const L1d_cache_miss_event_name = "perf_raw::rc860";
static auto const L2_cache_miss_event_name = "perf_raw::r0864";

bool init_papi() {
    if (int ret = PAPI_library_init(PAPI_VER_CURRENT); ret != PAPI_VER_CURRENT) {
        std::cerr << "Error initializing PAPI\n";
        return false;
    }
    if (int ret = PAPI_thread_init((unsigned long (*)())(pthread_self)); ret != PAPI_OK) {
        std::cerr << "Error initializing threads for PAPI\n";
        return false;
    }
    if (int ret = PAPI_query_named_event(L1d_cache_miss_event_name); ret != PAPI_OK) {
        std::cerr << "Cannot measure event '" << L1d_cache_miss_event_name << "'\n";
        return false;
    }
    if (int ret = PAPI_query_named_event(L2_cache_miss_event_name); ret != PAPI_OK) {
        std::cerr << "Cannot measure event '" << L2_cache_miss_event_name << "'\n";
        return false;
    }
    return true;
}

bool start_performance_counter(int& event_set) {
    if (int ret = PAPI_register_thread(); ret != PAPI_OK) {
        return false;
    }
    if (int ret = PAPI_create_eventset(&event_set); ret != PAPI_OK) {
        return false;
    }
    int event_code{};
    if (int ret = PAPI_event_name_to_code(L1d_cache_miss_event_name, &event_code); ret != PAPI_OK) {
        return false;
    }
    if (int ret = PAPI_add_event(event_set, event_code); ret != PAPI_OK) {
        return false;
    }
    if (int ret = PAPI_event_name_to_code(L2_cache_miss_event_name, &event_code); ret != PAPI_OK) {
        return false;
    }
    if (int ret = PAPI_add_event(event_set, event_code); ret != PAPI_OK) {
        return false;
    }
    if (int ret = PAPI_start(event_set); ret != PAPI_OK) {
        return false;
    }
    return true;
}
#endif

enum class WorkMode { Mixed, Split, Increment };
enum class KeyDistribution { Uniform, Ascending, Descending };

struct Settings {
    int num_threads = 4;
    std::size_t prefill_per_thread = 1 << 20;
    std::size_t elements_per_thread = 1 << 24;
    WorkMode work_mode = WorkMode::Mixed;
    int num_push_threads = 1;
    KeyDistribution key_distribution = KeyDistribution::Uniform;
#ifdef USE_UINT8
    key_type max_key = std::numeric_limits<key_type>::max();
#else
    key_type max_key = 1 << 30;
#endif
    int seed = 1;
#ifdef QUALITY
    std::filesystem::path histogram_file = "histogram.dat";
    std::filesystem::path log_file;
#endif
#ifdef WITH_PAPI
    bool enable_performance_counter = false;
#endif

    bool set_work_mode(char c) {
        switch (c) {
            case 'm':
                work_mode = WorkMode::Mixed;
                return true;
            case 's':
                work_mode = WorkMode::Split;
                return true;
            case 'i':
                work_mode = WorkMode::Increment;
                return true;
        }
        return false;
    }

    bool set_key_distribution(char c) {
        switch (c) {
            case 'u':
                key_distribution = KeyDistribution::Uniform;
                return true;
            case 'a':
                key_distribution = KeyDistribution::Ascending;
                return true;
            case 'd':
                key_distribution = KeyDistribution::Descending;
                return true;
        }
        return false;
    }

    [[nodiscard]] auto work_mode_str() const {
        switch (work_mode) {
            case WorkMode::Mixed:
                return "mixed";
            case WorkMode::Split:
                return "split";
            case WorkMode::Increment:
                return "increment";
        }
        return "unknown";
    }
    [[nodiscard]] auto key_distribution_str() const {
        switch (key_distribution) {
            case KeyDistribution::Uniform:
                return "uniform";
            case KeyDistribution::Ascending:
                return "ascending";
            case KeyDistribution::Descending:
                return "descending";
        }
        return "unknown";
    }

    [[nodiscard]] bool validate() const {
        if (num_threads <= 0) {
            return false;
        }
        if (max_key == 0) {
            return false;
        }
        if (work_mode == WorkMode::Split) {
            if ((num_push_threads < 0 || num_push_threads > num_threads) ||
                (num_push_threads == 0 && elements_per_thread > 0)) {
                return false;
            }
            if (num_push_threads == 0 && elements_per_thread > 0) {
                return false;
            }
        }
        return true;
    }
};

struct Result {
    std::pair<timepoint_type, timepoint_type> work_time{timepoint_type::max(), timepoint_type::min()};
    long long num_failed_pops{0};
    long long num_pops{0};
    long long num_pushes{0};
#ifdef WITH_PAPI
    long long l1d_cache_misses{0};
    long long l2_cache_misses{0};
#endif
#ifdef MQ_COUNT_STATS
    multiqueue::Counters counters;
#endif
    bool success{true};
};

struct GlobalData {
    Settings settings;
#ifdef QUALITY
    std::vector<quality::OperationLog> logs;
#endif
    Result result;
    std::mutex result_mutex;

    void update_result(const Result& r) {
        std::lock_guard<std::mutex> lock(result_mutex);
        result.work_time.first = std::min(result.work_time.first, r.work_time.first);
        result.work_time.second = std::max(result.work_time.second, r.work_time.second);
        result.num_failed_pops += r.num_failed_pops;
        result.num_pops += r.num_pops;
        result.num_pushes += r.num_pushes;
#ifdef WITH_PAPI
        result.l1d_cache_misses += r.l1d_cache_misses;
        result.l2_cache_misses += r.l2_cache_misses;
#endif
#ifdef MQ_COUNT_STATS
        result.counters += r.counters;
#endif
    }
};

struct ThreadData {
    thread_coordination::Context ctx;
    handle_type handle;
    Result result{};
#ifdef QUALITY
    quality::OperationLog logs;
#endif
};

void execute_mixed(ThreadData& thread_data, std::vector<key_type> const& keys) {
#ifdef QUALITY
    thread_data.push_log.reserve(thread_data.push_log.size() + keys.size());
    thread_data.pop_log.reserve(keys.size());
    auto value = quality::packed_value::pack(thread_data.ctx.get_id(), thread_data.push_log.size());
#endif
    thread_data.result.work_time =
        thread_data.ctx.execute_synchronized_blockwise(keys.size(), [&](auto start_index, auto count) {
            for (auto i = start_index; i < start_index + count; ++i) {
#ifdef QUALITY
                thread_data.handle.push({keys[i], value});
                auto push_tick = get_tick();
                thread_data.push_log.push_back({push_tick, keys[i], value});
                ++value;
#else
            thread_data.handle.push({keys[i], keys[i]});
#endif
                PriorityQueue::value_type retval;
#ifdef QUALITY
                auto pop_tick = get_tick();
#endif
                while (!thread_data.handle.try_pop(retval)) {
                    ++thread_data.result.num_failed_pops;
#ifdef QUALITY
                    pop_tick = get_tick();
#endif
                }
#ifdef QUALITY
                thread_data.pop_log.push_back({pop_tick, retval.first, retval.second});
#endif
            }
            thread_data.result.num_pops += count;
            thread_data.result.num_pushes += count;
        });
    thread_data.result.success = true;
}

void execute_split_push(ThreadData& thread_data, std::vector<key_type> const& keys) {
#ifdef QUALITY
    thread_data.push_log.reserve(thread_data.push_log.size() + keys.size());
    auto value = quality::packed_value::pack(ctx.get_id(), push_log.size());
#endif
    thread_data.result.work_time =
        thread_data.ctx.execute_synchronized_blockwise(keys.size(), [&](auto start_index, auto count) {
            for (auto i = start_index; i < start_index + count; ++i) {
#ifdef QUALITY
                thread_data.handle.push({keys[i], value});
                auto tick = get_tick();
                thread_data.push_log.push_back({tick, keys[i], value});
                ++value;
#else
                thread_data.handle.push({keys[i], keys[i]});
#endif
            }
            thread_data.result.num_pushes += count;
        });
    thread_data.result.success = true;
}

void execute_split_pop(ThreadData& thread_data, std::atomic_llong& pop_count) {
#ifdef QUALITY
    thread_data.pop_log.reserve(num_elements);
#endif
    thread_data.result.work_time = thread_data.ctx.execute_synchronized([&thread_data, &pop_count]() {
        do {
            long long num_pops = 0;
            PriorityQueue::value_type retval;
            while (thread_data.handle.try_pop(retval)) {
#ifdef QUALITY
                auto tick = get_tick();
                thread_data.pop_log.emplace_back(tick, static_cast<std::uint32_t>(retval.second));
#endif
                ++num_pops;
            }
            ++thread_data.result.num_failed_pops;
            if (num_pops == 0) {
                if (pop_count.load(std::memory_order_relaxed) == 0) {
                    break;
                }
            } else {
                thread_data.result.num_pops += num_pops;
                if (pop_count.fetch_sub(num_pops, std::memory_order_relaxed) - num_pops == 0) {
                    break;
                }
            }
        } while (true);
    });
    thread_data.result.success = true;
}

void execute_increment(ThreadData& thread_data, std::vector<key_type> const& keys) {
#ifdef QUALITY
    thread_data.push_log.reserve(push_log.size() + keys.size());
    thread_data.pop_log.reserve(keys.size());
    auto value = quality::packed_value::pack(ctx.get_id(), push_log.size());
#endif
    thread_data.result.work_time =
        thread_data.ctx.execute_synchronized_blockwise(keys.size(), [&](auto start_index, auto count) {
            for (auto i = start_index; i < start_index + count; ++i) {
                PriorityQueue::value_type retval;
                while (!thread_data.handle.try_pop(retval)) {
                    ++thread_data.result.num_failed_pops;
                }
#ifdef QUALITY
                auto tick = get_tick();
                thread_data.pop_log.emplace_back(tick, retval.second);
                handle.push({retval.first + keys[i], value});
                ++value;
                tick = get_tick();
                thread_data.push_log.push_back({tick, retval.first + keys[i]});
#else
            thread_data.handle.push({retval.first + keys[i], keys[i]});
#endif
            }
            thread_data.result.num_pops += count;
            thread_data.result.num_pushes += count;
        });
    thread_data.result.success = true;
}

void benchmark_thread(ThreadData& thread_data, std::vector<key_type>& keys, GlobalData& global_data) {
    std::seed_seq seed{global_data.settings.seed, thread_data.ctx.get_id()};
    std::default_random_engine rng(seed);

    if (thread_data.ctx.get_id() == 0) {
        std::clog << "Generating keys..." << std::flush;
    }

    auto first_key =
        keys.begin() + thread_data.ctx.get_id() * static_cast<std::ptrdiff_t>(global_data.settings.elements_per_thread);
    auto last_key = first_key + static_cast<std::ptrdiff_t>(global_data.settings.elements_per_thread);
    thread_data.ctx.execute_synchronized(
        [&rng, &global_data](auto first, auto last) {
            auto key_dist = std::uniform_int_distribution<key_type>(1, global_data.settings.max_key);
            std::generate(first, last, [&key_dist, &rng]() { return key_dist(rng); });
        },
        first_key, last_key);

    if (thread_data.ctx.get_id() == 0) {
        if (global_data.settings.key_distribution == KeyDistribution::Ascending) {
            std::sort(keys.begin(), keys.end());
        } else if (global_data.settings.key_distribution == KeyDistribution::Descending) {
            std::sort(keys.begin(), keys.end(), std::greater<>());
        }
        std::clog << "done\nPrefilling..." << std::flush;
    }

#ifdef QUALITY
    thread_data.push_log.reserve(global_data.settings.prefill_per_thread);
#endif

    if (global_data.settings.prefill_per_thread > 0) {
#ifdef QUALITY
        auto& push_log = shared_data.push_log[static_cast<std::size_t>(ctx.get_id())];
#endif
        thread_data.ctx.execute_synchronized([&thread_data, &rng, &global_data]() {
            auto key_dist = std::uniform_int_distribution<key_type>(1, global_data.settings.max_key);
            for (std::size_t i = 0; i < global_data.settings.prefill_per_thread; ++i) {
                auto key = key_dist(rng);
#ifdef QUALITY
                auto v = quality::packed_value::pack(ctx.get_id(), i);
                thread_data.handle.push({key, v});
                thread_data.push_log.push_back({{}, key});
#else
                thread_data.handle.push({key, key});
#endif
            }
        });
    }

#ifdef MQ_COUNT_STATS
    thread_data.handle.reset_counters();
#endif

    if (thread_data.ctx.get_id() == 0) {
        std::clog << "done\nRunning benchmark..." << std::flush;
    }

#ifdef WITH_PAPI
    bool papi_started = false;
    int event_set = PAPI_NULL;
    if (global_data.settings.enable_performance_counter) {
        papi_started = start_performance_counter(event_set);
        if (!papi_started) {
            thread_data.ctx.write(std::cerr) << "Failed to start counters\n";
            thread_data.result.success = false;
        }
    }
#endif

    switch (global_data.settings.work_mode) {
        case WorkMode::Mixed:
            execute_mixed(thread_data, keys);
            break;
        case WorkMode::Split:
            if (thread_data.ctx.get_id() < global_data.settings.num_push_threads) {
                execute_split_push(thread_data, keys);
            } else {
                static std::atomic_llong pop_count = static_cast<long long>(global_data.settings.prefill_per_thread +
                                                                            global_data.settings.elements_per_thread) *
                    global_data.settings.num_threads;
                execute_split_pop(thread_data, pop_count);
            }
            break;
        case WorkMode::Increment:
            execute_increment(thread_data, keys);
            break;
    }

#ifdef WITH_PAPI
    if (papi_started) {
        std::array<long long, 2> counter{};
        if (int ret = PAPI_stop(event_set, counter.data()); ret != PAPI_OK) {
            thread_data.ctx.write(std::cerr) << "Failed to stop counters\n";
            thread_data.result.success = false;
        } else {
            thread_data.result.l1d_cache_misses += counter[0];
            thread_data.result.l2_cache_misses += counter[1];
        }
    }
#endif
#ifdef MQ_COUNT_STATS
    thread_data.result.counters = thread_data.handle.get_counters();
#endif
    global_data.update_result(thread_data.result);
    if (thread_data.ctx.get_id() == 0) {
        std::clog << "done" << std::endl;
    }
}

void print_settings(Settings const& settings, PriorityQueueConfig const& pq_config) {
    std::cout << "Threads: " << settings.num_threads << '\n'
              << "Prefill per thread: " << settings.prefill_per_thread << '\n'
              << "Elements per thread: " << settings.elements_per_thread << '\n'
              << "Operation mode: " << settings.work_mode_str();
    if (settings.work_mode == WorkMode::Split) {
        std::cout << " (" << settings.num_push_threads << " push)";
    }
    std::cout << '\n';
    std::cout << "Key distribution: " << settings.key_distribution_str() << '\n'
              << "Max key: " << static_cast<unsigned long>(settings.max_key) << '\n'
              << "Seed: " << settings.seed << '\n';
    std::cout << "Priority queue: ";
    describe_pq(std::cout, pq_config);
    std::cout << '\n';
#ifdef USE_UINT8
    std::cout << "Data type: uint8_t\n";
#else
    std::cout << "Data type: unsigned long\n";
#endif
    std::cout << '\n' << '\n';
}

void print_global_data(GlobalData const& global_data) {
    std::cout << "Work time (s): " << std::setprecision(3)
              << std::chrono::duration<double>(global_data.result.work_time.second - global_data.result.work_time.first)
                     .count()
              << '\n';
    std::cout << "Failed pops: " << global_data.result.num_failed_pops << '\n';
#ifdef WITH_PAPI
    if (global_data.settings.enable_performance_counter) {
        std::cout << "L1d cache misses: " << global_data.result.l1d_cache_misses << '\n';
        std::cout << "L2 cache misses: " << global_data.result.l2_cache_misses << '\n';
    }
#endif
}

void print_csv(Settings const& settings, GlobalData const& global_data, std::filesystem::path const& file) {
    std::ofstream out(file);
    if (!out) {
        throw std::runtime_error("Failed to open file: " + file.string());
    }
    out << "threads,prefill,elements,work-mode,push-threads,key-distribution,max-key,seed,work-time,"
           "failed-pops,l1d-cache-misses,l2-cache-misses,num-resets,use-counts,success\n";
    out << settings.num_threads << ',' << settings.prefill_per_thread << ',' << settings.elements_per_thread << ','
        << settings.work_mode_str() << ',' << settings.num_push_threads << ',' << settings.key_distribution_str() << ','
        << static_cast<unsigned long>(settings.max_key) << ',' << settings.seed << ',' << std::fixed
        << std::setprecision(3)
        << std::chrono::duration<double>(global_data.result.work_time.second - global_data.result.work_time.first)
               .count()
        << ',' << global_data.result.num_failed_pops;
#ifdef WITH_PAPI
    if (global_data.settings.enable_performance_counter) {
        out << ',' << global_data.result.l1d_cache_misses << ',' << global_data.result.l2_cache_misses;
    } else {
        out << ",n/a,n/a";
    }
#else
    out << ",n/a,n/a";
#endif
#ifdef MQ_COUNT_STATS
    out << ",n/a,n/a";
#else
    out << ",n/a,n/a";
#endif
    out << ',' << global_data.result.success << std::endl;
}

void print_header() {
    std::clog << "Built on " << __DATE__ << ' ' << __TIME__ << " with:\n";
#ifdef NDEBUG
    std::clog << "  NDEBUG defined\n";
#else
    std::clog << "  NDEBUG not defined\n";
#endif
#ifdef __clang__
    std::clog << "  Clang " << __clang_version__ << '\n';
#elif defined(__GNUC__)
    std::clog << "  GCC " << __VERSION__ << '\n';
#else
    std::clog << "  Unknown compiler\n";
#endif
#ifdef WITH_PAPI
    std::clog << "  PAPI " << PAPI_VER_CURRENT << '\n';
#else
    std::clog << "  PAPI disabled\n";
#endif
#ifdef MQ_COUNT_STATS
    std::clog << "  MQ_COUNT_STATS enabled\n";
#else
    std::clog << "  MQ_COUNT_STATS disabled\n";
#endif
}

int main(int argc, char* argv[]) {
    print_header();

    std::clog << "\nCommand line:";
    for (int i = 0; i < argc; ++i) {
        std::clog << ' ' << argv[i];
    }
    std::clog << '\n' << '\n';

    GlobalData global_data;
    cxxopts::Options options(argv[0]);
    // clang-format off
    options.add_options()
        ("h,help", "Print this help")
        ("j,threads", "The number of threads", cxxopts::value<int>(global_data.settings.num_threads), "NUMBER")
        ("p,prefill", "The prefill per thread", cxxopts::value<std::size_t>(global_data.settings.prefill_per_thread), "NUMBER")
        ("n,keys", "The number of keys per thread", cxxopts::value<std::size_t>(global_data.settings.elements_per_thread), "NUMBER")
        ("w,work-mode", "Specify the work mode ([m]ixed, [s]plit, [i]ncrement)", cxxopts::value<char>(), "STRING")
        ("i,push-threads", "The number of pushing threads in split mode", cxxopts::value<int>(global_data.settings.num_push_threads), "NUMBER")
        ("e,element-distribution", "Specify the element distribution ([u]niform, [a]scending, [d]escending)", cxxopts::value<char>(), "STRING")
        ("m,max", "Specify the max key", cxxopts::value<key_type>(global_data.settings.max_key), "NUMBER")
        ("s,seed", "Specify the initial seed", cxxopts::value<int>(global_data.settings.seed), "NUMBER")
#ifdef QUALITY
        ("q,histogram-file", "Path to write the histogram to", cxxopts::value<std::filesystem::path>(settings.histogram_file), "PATH")
        ("x,log-file", "Path to write the operation log to", cxxopts::value<std::filesystem::path>(settings.log_file), "PATH")
#endif
#ifdef WITH_PAPI
        ("r,pc", "Enable performance counter", cxxopts::value<>(global_data.settings.enable_performance_counter))
#endif
        ;
    // clang-format on
    add_pq_options(options);

    PriorityQueueConfig pq_config;

    {
        cxxopts::ParseResult result;
        try {
            result = options.parse(argc, argv);
        } catch (cxxopts::OptionParseException const& e) {
            std::cerr << "Error parsing arguments: " << e.what() << '\n';
            std::cerr << options.help() << std::endl;
            return 1;
        }
        if (result.count("help") > 0) {
            std::cerr << options.help() << std::endl;
            return 0;
        }
        if (result.count("work-mode") > 0) {
            auto work_mode = result["work-mode"].as<char>();
            if (!global_data.settings.set_work_mode(work_mode)) {
                std::cerr << "Invalid work mode: " << result["work-mode"].as<char>() << '\n';
                std::cerr << options.help() << std::endl;
                return 1;
            }
        }
        if (result.count("element-distribution") > 0) {
            auto key_distribution = result["element-distribution"].as<char>();
            if (!global_data.settings.set_key_distribution(key_distribution)) {
                std::cerr << "Invalid element distribution: " << result["element-distribution"].as<char>() << '\n';
                std::cerr << options.help() << std::endl;
                return 1;
            }
        }
        pq_config = get_pq_options(result);
    }

    print_settings(global_data.settings, pq_config);

    if (!global_data.settings.validate()) {
        std::cerr << "Error: Invalid settings\n";
        std::cerr << options.help() << std::endl;
        return 1;
    }

    auto pq = create_pq<PriorityQueue>(
        global_data.settings.num_threads,
        global_data.settings.prefill_per_thread * static_cast<std::size_t>(global_data.settings.num_threads),
        pq_config);
    std::vector<key_type> keys(static_cast<std::size_t>(global_data.settings.num_threads) *
                               global_data.settings.elements_per_thread);
#ifdef WITH_PAPI
    if (global_data.settings.enable_performance_counter) {
        if (!init_papi()) {
            global_data.result.success = false;
        }
    }
#endif

    thread_coordination::TaskHandle task_handle(global_data.settings.num_threads, [&keys, &global_data, &pq](auto ctx) {
        ThreadData thread_data{ctx, pq.get_handle()};
        benchmark_thread(thread_data, keys, global_data);
    });
    task_handle.wait();

#ifdef QUALITY
    if (!settings.log_file.empty()) {
        std::clog << "Writing logs..." << std::flush;
        try {
            quality::write_logs(shared_data.push_log, shared_data.pop_log, settings.log_file);
            std::clog << "done" << std::endl;
        } catch (std::exception const& e) {
            std::clog << "failed\n" << std::endl;
            std::cerr << "Error: " << e.what() << std::endl;
            shared_data.success = false;
        }
    }
    std::clog << "Checking logs..." << std::flush;
    try {
        quality::fix_logs(shared_data.push_log, shared_data.pop_log);
        std::clog << "done" << std::endl;
    } catch (std::runtime_error const& e) {
        std::clog << "failed" << std::endl;
        std::cerr << "Error: " << e.what() << std::endl;
        shared_data.success = false;
    }
    if (shared_data.success && !settings.histogram_file.empty()) {
        try {
            quality::write_histogram(shared_data.push_log, shared_data.pop_log, settings.histogram_file);
            std::clog << "done" << std::endl;
        } catch (std::runtime_error const& e) {
            std::clog << "failed" << std::endl;
            std::err << "Error: " << e.what() << std::endl;
            shared_data.success = false;
        }
    }
#endif

    std::clog << '\n';
    print_global_data(global_data);
}
