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
using handle_type = PriorityQueue::Handle;

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

struct Settings {
    enum class WorkMode { Mixed, Split, Increment };
    enum class ElementDistribution { Uniform, Ascending, Descending };

    int num_threads = 4;
    std::size_t prefill_per_thread = 1 << 20;
    std::size_t elements_per_thread = 1 << 24;
    WorkMode work_mode = WorkMode::Mixed;
    int num_push_threads = 1;
    ElementDistribution element_distribution = ElementDistribution::Uniform;
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

    bool set_element_distribution(char c) {
        switch (c) {
            case 'u':
                element_distribution = ElementDistribution::Uniform;
                return true;
            case 'a':
                element_distribution = ElementDistribution::Ascending;
                return true;
            case 'd':
                element_distribution = ElementDistribution::Descending;
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
    [[nodiscard]] auto element_distribution_str() const {
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

Settings& get_settings() {
    static Settings settings;
    return settings;
}

struct Result {
    timepoint_type start_time{timepoint_type::max()};
    timepoint_type end_time{timepoint_type::min()};
    long long num_failed_pops{0};
    long long num_pops{0};
#ifdef WITH_PAPI
    long long l1d_cache_misses{0};
    long long l2_cache_misses{0};
#endif
#ifdef MQ_COUNT_STATS
    long long num_locking_failed{0};
    long long num_resets{0};
    long long use_counts{0};
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

    void update_result(const Result& result) {
        auto old_start = start_time.load(std::memory_order_relaxed);
        while (work_time.first < old_start &&
               !start_time.compare_exchange_weak(old_start, work_time.first, std::memory_order_relaxed)) {
        }
        auto old_end = start_time.load(std::memory_order_relaxed);
        while (work_time.second > old_end &&
               !end_time.compare_exchange_weak(old_end, work_time.second, std::memory_order_relaxed)) {
        }
    }
};

template <typename RNG>
void generate_workload(std::vector<key_type>& keys, int id, RNG& rng) {
    auto start_index = static_cast<std::size_t>(id) * get_settings().elements_per_thread;
    switch (get_settings().element_distribution) {
        case Settings::ElementDistribution::Uniform: {
            auto key_dist = std::uniform_int_distribution<key_type>(1, get_settings().max_key);
            for (std::size_t i = start_index; i < start_index + get_settings().elements_per_thread; ++i) {
                keys[i] = key_dist(rng);
            }
        } break;
        case Settings::ElementDistribution::Ascending: {
            auto m = get_settings().max_key;
            for (std::size_t i = start_index; i < start_index + get_settings().elements_per_thread; ++i) {
                keys[i] = static_cast<key_type>((i * m) / keys.size() + 1);
            }
        } break;
        case Settings::ElementDistribution::Descending: {
            auto m = get_settings().max_key;
            for (std::size_t i = start_index; i < start_index + get_settings().elements_per_thread; ++i) {
                keys[i] = static_cast<key_type>(((keys.size() - i - 1) * m) / keys.size() + 1);
            }
        } break;
    }
}

void execute_mixed(thread_coordination::Context const& ctx, handle_type& handle, std::vector<key_type> const& keys,
                   SharedData& shared_data) {
    long long num_failed_pops = 0;
#ifdef QUALITY
    auto& push_log = shared_data.push_log[static_cast<std::size_t>(ctx.get_id())];
    push_log.reserve(push_log.size() + keys.size());
    auto& pop_log = shared_data.pop_log[static_cast<std::size_t>(ctx.get_id())];
    pop_log.reserve(keys.size());
    auto value = quality::packed_value::pack(ctx.get_id(), push_log.size());
#endif
    auto work_time = ctx.execute_synchronized_blockwise(keys.size(), [&](auto start_index, auto count) {
        for (auto i = start_index; i < start_index + count; ++i) {
#ifdef QUALITY
            handle.push({keys[i], value});
            auto push_tick = get_tick();
            push_log.push_back({push_tick, keys[i], value});
            ++value;
#else
            handle.push({keys[i], keys[i]});
#endif
            PriorityQueue::value_type retval;
#ifdef QUALITY
            auto pop_tick = get_tick();
#endif
            while (!handle.try_pop(retval)) {
                ++num_failed_pops;
#ifdef QUALITY
                pop_tick = get_tick();
#endif
            }
#ifdef QUALITY
            pop_log.push_back({pop_tick, retval.first, retval.second});
#endif
        }
    });
    shared_data.num_failed_pops += num_failed_pops;
    shared_data.update_work_time(work_time);
}

void execute_split_push(thread_coordination::Context const& ctx, handle_type& handle, std::vector<key_type> const& keys,
                        SharedData& shared_data) {
#ifdef QUALITY
    auto& push_log = shared_data.push_log[static_cast<std::size_t>(ctx.get_id())];
    push_log.reserve(push_log.size() + keys.size());
    auto index = push_log.size();
    auto value = quality::packed_value::pack(ctx.get_id(), push_log.size());
#endif
    auto work_time = ctx.execute_synchronized_blockwise(keys.size(), [&](auto start_index, auto count) {
        for (auto i = start_index; i < start_index + count; ++i) {
#ifdef QUALITY
            handle.push({keys[i], value});
            auto tick = get_tick();
            push_log.push_back({tick, keys[i], value});
            ++value;
#else
            handle.push({keys[i], keys[i]});
#endif
        }
    });
    shared_data.update_work_time(work_time);
}

void execute_split_pop(thread_coordination::Context const& ctx, handle_type& handle, SharedData& shared_data,
                       std::size_t num_elements) {
    long long num_failed_pops = 0;
#ifdef QUALITY
    auto& pop_log = shared_data.pop_log[static_cast<std::size_t>(ctx.get_id())];
    pop_log.reserve(num_elements);
#endif
    auto work_time = ctx.execute_synchronized([&, num_elements]() {
        do {
            long long num_pops = 0;
            PriorityQueue::value_type retval;
            while (handle.try_pop(retval)) {
#ifdef QUALITY
                auto tick = get_tick();
                pop_log.emplace_back(tick, static_cast<std::uint32_t>(retval.second));
#endif
                ++num_pops;
            }
            ++num_failed_pops;
            if (num_pops == 0) {
                if (shared_data.num_pops.load(std::memory_order_relaxed) >= static_cast<long long>(num_elements)) {
                    break;
                }
            } else {
                if (shared_data.num_pops.fetch_add(num_pops, std::memory_order_relaxed) + num_pops >=
                    static_cast<long long>(num_elements)) {
                    break;
                }
            }
        } while (true);
    });
    assert(shared_data.num_pops.load(std::memory_order_relaxed) == static_cast<long long>(num_elements));
    shared_data.update_work_time(work_time);
    shared_data.num_failed_pops += num_failed_pops;
}

void execute_increment(thread_coordination::Context const& ctx, handle_type& handle, std::vector<key_type> const& keys,
                       SharedData& shared_data) {
    long long num_failed_pops = 0;
#ifdef QUALITY
    auto& push_log = shared_data.push_log[static_cast<std::size_t>(ctx.get_id())];
    push_log.reserve(push_log.size() + keys.size());
    auto& pop_log = shared_data.pop_log[static_cast<std::size_t>(ctx.get_id())];
    pop_log.reserve(keys.size());
    auto index = push_log.size();
#endif
    auto work_time = ctx.execute_synchronized_blockwise(keys.size(), [&](auto start_index, auto count) {
        for (auto i = start_index; i < start_index + count; ++i) {
            PriorityQueue::value_type retval;
            while (!handle.try_pop(retval)) {
                ++num_failed_pops;
            }
#ifdef QUALITY
            auto tick = get_tick();
            pop_log.emplace_back(tick, retval.second);
            auto v = quality::packed_value::pack(ctx.get_id(), index);
            handle.push({retval.first + keys[i], v});
            ++index;
            tick = get_tick();
            push_log.push_back({tick, retval.first + keys[i]});
#else
            handle.push({retval.first + keys[i], keys[i]});
#endif
        }
    });
    shared_data.num_failed_pops += num_failed_pops;
    shared_data.update_work_time(work_time);
}

void benchmark_thread(thread_coordination::Context ctx, PriorityQueue& pq, std::vector<key_type>& keys,
                      SharedData& shared_data) {
    std::seed_seq seed{get_settings().seed, ctx.get_id()};
    std::default_random_engine rng(seed);

    if (ctx.get_id() == 0) {
        std::clog << "Generating keys..." << std::flush;
    }
    ctx.execute_synchronized([&, id = ctx.get_id()]() { generate_workload(keys, id, rng); });

    if (ctx.get_id() == 0) {
        std::clog << "done\nPrefilling..." << std::flush;
    }

#ifdef QUALITY
    shared_data.push_log[static_cast<std::size_t>(ctx.get_id())].reserve(settings.prefill_per_thread);
#endif

    handle_type handle = pq.get_handle(ctx.get_id());

    if (get_settings().prefill_per_thread > 0) {
#ifdef QUALITY
        auto& push_log = shared_data.push_log[static_cast<std::size_t>(ctx.get_id())];
#endif
        auto key_dist = std::uniform_int_distribution<key_type>(1, get_settings().max_key);
        ctx.execute_synchronized([&, n = get_settings().prefill_per_thread]() {
            for (std::size_t i = 0; i < n; ++i) {
                auto key = key_dist(rng);
#ifdef QUALITY
                auto v = quality::packed_value::pack(ctx.get_id(), i);
                handle.push({key, v});
                push_log.push_back({{}, key});
#else
                handle.push({key, key});
#endif
            }
        });
    }

#ifdef MQ_COUNT_STATS
    handle.stats.reset();
#endif

    if (ctx.get_id() == 0) {
        std::clog << "done\nRunning benchmark..." << std::flush;
    }

#ifdef WITH_PAPI
    bool papi_started = false;
    int event_set = PAPI_NULL;
    if (settings.enable_performance_counter) {
        papi_started = start_performance_counter(event_set);
        if (!papi_started) {
            ctx.write(std::cerr) << "Failed to start counters\n";
            shared_data.success = false;
        }
    }
#endif

    switch (get_settings().work_mode) {
        case Settings::WorkMode::Mixed:
            execute_mixed(ctx, handle, keys, shared_data);
            break;
        case Settings::WorkMode::Split:
            if (ctx.get_id() < get_settings().num_push_threads) {
                execute_split_push(ctx, handle, keys, shared_data);
            } else {
                execute_split_pop(ctx, handle, shared_data,
                                  (get_settings().prefill_per_thread + get_settings().elements_per_thread) *
                                      static_cast<std::size_t>(get_settings().num_threads));
            }
            break;
        case Settings::WorkMode::Increment:
            execute_increment(ctx, handle, keys, shared_data);
            break;
    }

#ifdef WITH_PAPI
    if (papi_started) {
        std::array<long long, 2> counter{};
        if (int ret = PAPI_stop(event_set, counter.data()); ret != PAPI_OK) {
            ctx.write(std::cerr) << "Failed to stop counters\n";
            shared_data.success = false;
        } else {
            shared_data.l1d_cache_misses += counter[0];
            shared_data.l2_cache_misses += counter[1];
        }
    }
#endif

#ifdef MQ_COUNT_STATS
    shared_data.num_locking_failed += handle.stats.num_locking_failed;
    shared_data.num_resets += handle.stats.num_resets;
    shared_data.use_counts += handle.stats.use_counts;
#endif

    if (ctx.get_id() == 0) {
        std::clog << "done" << std::endl;
    }
}

void print_settings(Settings const& settings, PriorityQueueConfig const& pq_config) {
    std::cout << "Threads: " << settings.num_threads << '\n'
              << "Prefill per thread: " << settings.prefill_per_thread << '\n'
              << "Elements per thread: " << settings.elements_per_thread << '\n'
              << "Operation mode: " << settings.work_mode_str();
    if (settings.work_mode == Settings::WorkMode::Split) {
        std::cout << " (" << settings.num_push_threads << " push)";
    }
    std::cout << '\n';
    std::cout << "Element distribution: " << settings.element_distribution_str() << '\n'
              << "Max key: " << static_cast<unsigned long>(settings.max_key) << '\n'
              << "Seed: " << settings.seed << '\n';
    std::cout << "Priority queue: " << pq_name << ' ';
    describe_pq(std::cout, pq_config);
#ifdef USE_UINT8
    std::cout << "Data type: uint8_t\n";
#else
    std::cout << "Data type: unsigned long\n";
#endif
    std::cout << '\n' << '\n';
}

void print_shared_data(SharedData const& shared_data) {
    std::cout << "Work time (s): " << std::setprecision(3)
              << std::chrono::duration<double>(shared_data.end_time.load() - shared_data.start_time.load()).count()
              << '\n';
    std::cout << "Failed pops: " << shared_data.num_failed_pops << '\n';
#ifdef WITH_PAPI
    if (settings.enable_performance_counter) {
        std::cout << "L1d cache misses: " << shared_data.l1d_cache_misses << '\n';
        std::cout << "L2 cache misses: " << shared_data.l2_cache_misses << '\n';
    }
#endif
#ifdef MQ_COUNT_STATS
    std::cout << "Failed locks per operation: "
              << static_cast<double>(shared_data.num_locking_failed) /
            static_cast<double>(static_cast<std::size_t>(settings.num_threads) * settings.elements_per_thread)
              << '\n';
    std::cout << "Average queue use count: "
              << static_cast<double>(shared_data.use_counts) / static_cast<double>(shared_data.num_resets) << '\n';
#endif
}

void print_csv(Settings const& settings, SharedData const& shared_data, std::filesystem::path const& file) {
    std::ofstream out(file);
    if (!out) {
        throw std::runtime_error("Failed to open file: " + file.string());
    }
    out << "threads,prefill,elements,work-mode,push-threads,element-distribution,max-key,seed,work-time,"
           "failed-pops,l1d-cache-misses,l2-cache-misses,num-resets,use-counts,success\n";
    out << settings.num_threads << ',' << settings.prefill_per_thread << ',' << settings.elements_per_thread << ','
        << settings.work_mode_str() << ',' << settings.num_push_threads << ',' << settings.element_distribution_str()
        << ',' << static_cast<unsigned long>(settings.max_key) << ',' << settings.seed << ',' << std::fixed
        << std::setprecision(3)
        << std::chrono::duration<double>(shared_data.end_time.load() - shared_data.start_time.load()).count() << ','
        << shared_data.num_failed_pops;
#ifdef WITH_PAPI
    if (settings.enable_performance_counter) {
        out << ',' << shared_data.l1d_cache_misses << ',' << shared_data.l2_cache_misses;
    } else {
        out << ",n/a,n/a";
    }
#else
    out << ",n/a,n/a";
#endif
#ifdef MQ_COUNT_STATS
    out << ',' << shared_data.num_resets << ',' << shared_data.use_counts;
#else
    out << ",n/a,n/a";
#endif
    out << ',' << shared_data.success.load() << std::endl;
}

bool run_benchmark(PriorityQueueConfig const& pq_config) {
    auto pq = create_pq<PriorityQueue>(
        get_settings().num_threads,
        get_settings().prefill_per_thread * static_cast<std::size_t>(get_settings().num_threads), pq_config);
    std::vector<key_type> keys(static_cast<std::size_t>(get_settings().num_threads) *
                               get_settings().elements_per_thread);
    SharedData shared_data;
#ifdef QUALITY
    shared_data.pop_log.resize(static_cast<std::size_t>(settings.num_threads));
    shared_data.push_log.resize(static_cast<std::size_t>(settings.num_threads));
#endif

#ifdef WITH_PAPI
    if (settings.enable_performance_counter) {
        if (!init_papi()) {
            shared_data.success = false;
        }
    }
#endif

    thread_coordination::TaskHandle task_handle(settings.num_threads, benchmark_thread, settings, std::ref(pq),
                                                std::ref(keys), std::ref(shared_data));
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
    print_shared_data(shared_data);
    return shared_data.success.load();
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

    Settings settings{};
    cxxopts::Options options(argv[0]);
    // clang-format off
    options.add_options()
        ("h,help", "Print this help")
        ("j,threads", "The number of threads", cxxopts::value<int>(settings.num_threads), "NUMBER")
        ("p,prefill", "The prefill per thread", cxxopts::value<std::size_t>(settings.prefill_per_thread), "NUMBER")
        ("n,keys", "The number of keys per thread", cxxopts::value<std::size_t>(settings.elements_per_thread), "NUMBER")
        ("w,work-mode", "Specify the work mode ([m]ixed, [s]plit, [i]ncrement)", cxxopts::value<char>(), "STRING")
        ("i,push-threads", "The number of pushing threads in split mode", cxxopts::value<int>(settings.num_push_threads), "NUMBER")
        ("e,element-distribution", "Specify the element distribution ([u]niform, [a]scending, [d]escending)", cxxopts::value<char>(), "STRING")
        ("m,max", "Specify the max key", cxxopts::value<key_type>(settings.max_key), "NUMBER")
        ("s,seed", "Specify the initial seed", cxxopts::value<int>(settings.seed), "NUMBER")
#ifdef QUALITY
        ("q,histogram-file", "Path to write the histogram to", cxxopts::value<std::filesystem::path>(settings.histogram_file), "PATH")
        ("x,log-file", "Path to write the operation log to", cxxopts::value<std::filesystem::path>(settings.log_file), "PATH")
#endif
#ifdef WITH_PAPI
        ("r,pc", "Enable performance counter", cxxopts::value<>(settings.enable_performance_counter))
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
            if (!settings.set_work_mode(work_mode)) {
                std::cerr << "Invalid work mode: " << result["work-mode"].as<char>() << '\n';
                std::cerr << options.help() << std::endl;
                return 1;
            }
        }
        if (result.count("element-distribution") > 0) {
            auto element_distribution = result["element-distribution"].as<char>();
            if (!settings.set_element_distribution(element_distribution)) {
                std::cerr << "Invalid element distribution: " << result["element-distribution"].as<char>() << '\n';
                std::cerr << options.help() << std::endl;
                return 1;
            }
        }
        pq_config = get_pq_options(result);
    }

    print_settings(settings, pq_config);

    if (!settings.validate()) {
        std::cerr << "Error: Invalid settings\n";
        std::cerr << options.help() << std::endl;
        return 1;
    }

    bool success = run_benchmark(settings, pq_config);

    return success ? 0 : 1;
}
