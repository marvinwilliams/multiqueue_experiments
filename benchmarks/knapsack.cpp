#include "knapsack_instance.hpp"
#include "priority_queue_factory.hpp"
#include "termination_detection.hpp"
#include "thread_coordination.hpp"

#include "cxxopts.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <new>
#include <numeric>
#include <random>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

#ifdef MQ_HAS_MAX

using PriorityQueue = DefaultMaxPriorityQueue;

void push(PriorityQueue& pq, PriorityQueue::value_type const& v) {
    pq.push(v);
}

bool try_pop(PriorityQueue& pq, PriorityQueue::value_type& v) {
    return pq.try_pop(v);
}

#else

using PriorityQueue = DefaultMinPriorityQueue;

static constexpr unsigned long max_key = 1 << 28;

void push(PriorityQueue::Handle& h, PriorityQueue::value_type v) {
    h.push({max_key - v.first, v.second});
}

bool try_pop(PriorityQueue::Handle& h, PriorityQueue::value_type& v) {
    bool success = h.try_pop(v);
    if (!success) {
        return false;
    }
    v.first = max_key - v.first;
    return true;
}

#endif

using handle_type = PriorityQueue::Handle;

using thread_coordination::timepoint_type;

static_assert(sizeof(unsigned long) >= sizeof(std::uint64_t), "64bit unsigned long required");
using payload_type = unsigned long;

static constexpr std::uint8_t bits_for_index = 16;
static constexpr std::uint8_t bits_for_free_capacity = 24;
static constexpr std::uint8_t bits_for_value = 24;
static_assert(bits_for_index + bits_for_free_capacity + bits_for_value <= std::numeric_limits<payload_type>::digits,
              "Too many bits required for payload");

constexpr payload_type to_payload(std::size_t index, unsigned long free_capacity, unsigned long value) noexcept {
    assert(index < (1UL << bits_for_index));
    assert(free_capacity < (1UL << bits_for_free_capacity));
    assert(value < (1UL << bits_for_value));
    payload_type payload = value;
    payload <<= bits_for_free_capacity;
    payload |= free_capacity;
    payload <<= bits_for_index;
    payload |= index;
    return payload;
}

struct Settings {
    int num_threads = 4;
    std::filesystem::path knapsack_file;
    int seed = 1;
    unsigned long solution = 0;
};

struct ThreadData {
    handle_type pq_handle;
    long long pushed_nodes{0};
    long long ignored_nodes{0};
    long long processed_nodes{0};
};

struct SharedData {
    alignas(L1_CACHE_LINESIZE) std::atomic_ulong best_value{0};
    std::atomic<timepoint_type> start_time{timepoint_type::max()};
    std::atomic<timepoint_type> end_time{timepoint_type::min()};
    std::atomic_llong pushed_nodes{0};
    std::atomic_llong ignored_nodes{0};
    std::atomic_llong processed_nodes{0};
    termination_detection::Data termination_detection_data;

    void update_work_time(std::pair<timepoint_type, timepoint_type> const& work_time) {
        auto old_start = start_time.load(std::memory_order_relaxed);
        while (work_time.first < old_start &&
               !start_time.compare_exchange_weak(old_start, work_time.first, std::memory_order_relaxed)) {
        }
        auto old_end = start_time.load(std::memory_order_relaxed);
        while (work_time.second > old_end &&
               !end_time.compare_exchange_weak(old_end, work_time.second, std::memory_order_relaxed)) {
        }
    }

    void update_counters(ThreadData const& thread_data) {
        pushed_nodes.fetch_add(thread_data.pushed_nodes, std::memory_order_relaxed);
        ignored_nodes.fetch_add(thread_data.ignored_nodes, std::memory_order_relaxed);
        processed_nodes.fetch_add(thread_data.processed_nodes, std::memory_order_relaxed);
    }
};

unsigned long lower_bound(KnapsackInstance const& instance, unsigned long capacity, std::size_t index) noexcept {
    unsigned long value{0};
    while (index < instance.items.size() && instance.items[index].weight <= capacity) {
        capacity -= instance.items[index].weight;
        value += instance.items[index].value;
        ++index;
    }
    return value;
}

unsigned long upper_bound(KnapsackInstance const& instance, unsigned long capacity, std::size_t index) noexcept {
    assert(index <= instance.items.size());
    unsigned long value_offset = instance.prefix_sum[index].value;
    unsigned long target_capacity = instance.prefix_sum[index].weight + capacity;
    while (index != instance.prefix_sum.size()) {
        if (instance.prefix_sum[index].weight > target_capacity) {
            double fraction = static_cast<double>(target_capacity - instance.prefix_sum[index - 1].weight) /
                static_cast<double>(instance.items[index - 1].weight);
            return (instance.prefix_sum[index - 1].value - value_offset) +
                static_cast<unsigned long>(static_cast<double>(instance.items[index - 1].value) * fraction);
        }
        ++index;
    }
    // All items fit
    return instance.prefix_sum[index - 1].value - value_offset;
}

void process_node(PriorityQueue::value_type const& node, SharedData& shared_data, ThreadData& data,
                  KnapsackInstance const& instance) {
    auto best_value = shared_data.best_value.load(std::memory_order_relaxed);
    if (node.first <= best_value) {
        // The upper bound of this node is worse than the currently best value
        ++data.ignored_nodes;
        return;
    }
    ++data.processed_nodes;
    auto e = node.second;
    std::size_t index = e & ((1UL << bits_for_index) - 1);
    e >>= bits_for_index;
    unsigned long free_capacity = e & ((1UL << bits_for_free_capacity) - 1);
    assert(free_capacity <= instance.capacity);
    e >>= bits_for_free_capacity;
    unsigned long value = e & ((1UL << bits_for_value) - 1);
    if (index == instance.items.size()) {
        while (value > best_value &&
               !shared_data.best_value.compare_exchange_weak(best_value, value, std::memory_order_relaxed)) {
        }
        return;
    }
    /* std::cout << "popped " << index << ' ' << hint << ' ' << free_capacity << ' ' << value << std::endl; */
    if (instance.items[index].weight <= free_capacity) {
        auto new_value = value + instance.items[index].value;
        auto new_capacity = free_capacity - instance.items[index].weight;
        /* std::cout << "pushed " << index + 1 << ' ' << hint << ' ' << new_capacity << ' ' << new_value <<
         * std::endl;
         */
        auto payload = to_payload(index + 1, new_capacity, new_value);
        push(data.pq_handle, {node.first, payload});
        ++data.pushed_nodes;
    }
    auto upper_bound_without_next = value + upper_bound(instance, free_capacity, index + 1);
    if (upper_bound_without_next <= best_value) {
        return;
    }
    push(data.pq_handle, {upper_bound_without_next, to_payload(index + 1, free_capacity, value)});
    ++data.pushed_nodes;
}

void main_loop(int num_threads, SharedData& shared_data, ThreadData& data, KnapsackInstance const& instance) {
    PriorityQueue::value_type retval;
    while (termination_detection::try_do(num_threads, shared_data.termination_detection_data,
                                         [&]() { return try_pop(data.pq_handle, retval); })) {
        process_node(retval, shared_data, data, instance);
    }
}

void benchmark_thread(thread_coordination::Context ctx, PriorityQueue& pq, SharedData& shared_data,
                      KnapsackInstance const& instance) {
    ThreadData data{pq.get_handle(ctx.get_id())};
    if (ctx.get_id() == 0) {
        auto ub = upper_bound(instance, instance.capacity, 0);
        auto lb = lower_bound(instance, instance.capacity, 0);
        shared_data.best_value.store(lb, std::memory_order_relaxed);
        /* std::cout << "pushed " << 0 << ' ' << hint << ' ' << instance.capacity << ' ' << 0 << std::endl; */
        push(data.pq_handle, {ub, to_payload(0, instance.capacity, 0)});
        ++data.pushed_nodes;
        std::clog << "Solving knapsack instance..." << std::flush;
    }
    auto work_time = ctx.execute_synchronized([&]() { main_loop(ctx.get_num_threads(), shared_data, data, instance); });
    shared_data.update_work_time(work_time);
    shared_data.update_counters(data);
}

void print_settings(Settings const& settings) {
    std::clog << "Threads: " << settings.num_threads << '\n'
              << "Instance: " << settings.knapsack_file.string() << '\n'
              << "Seed: " << settings.seed;
    std::clog << "\n\n";
}

void print_shared_data(Settings const& settings, SharedData const& shared_data, KnapsackInstance const& instance,
                       bool valid) {
    auto time = std::chrono::duration<double>(shared_data.end_time.load() - shared_data.start_time.load()).count();
    std::clog << "Time (s): " << std::setprecision(3) << time << '\n';
    std::clog << "Value: " << shared_data.best_value.load() << '\n';
    std::clog << "Processed nodes: " << shared_data.processed_nodes.load() << '\n';
    std::clog << "Ignored nodes: " << shared_data.ignored_nodes.load() << '\n';
    std::clog << "Total nodes: " << shared_data.processed_nodes.load() + shared_data.ignored_nodes.load() << '\n';

    std::cout << "instance,items,threads,seed,work_time,processed_nodes,ignored_nodes,value,valid\n";
    std::cout << settings.knapsack_file.string() << ',' << instance.items.size() << ',' << settings.num_threads << ','
              << settings.seed << ',' << time << ',' << shared_data.processed_nodes.load() << ','
              << shared_data.ignored_nodes << ',' << shared_data.best_value.load() << ',' << valid << '\n';
}

bool run_benchmark(Settings const& settings, PriorityQueueConfig const& pq_config) {
    std::clog << "Reading instance..." << std::flush;
    KnapsackInstance instance;
    try {
        instance.from_file(settings.knapsack_file);
    } catch (std::runtime_error const& e) {
        std::clog << "failed: " << e.what() << std::endl;
        return false;
    }

    if (instance.items.size() + 1 >= (1 << bits_for_index) || instance.capacity >= (1 << bits_for_free_capacity) ||
        settings.solution >= (1 << bits_for_value)) {
        std::clog << "failed: Instance cannot be represented\n";
        return false;
    }
    std::clog << "done\n";

    auto pq = create_pq<PriorityQueue>(settings.num_threads, 1UL << 20, pq_config);
    SharedData shared_data;
    thread_coordination::TaskHandle task_handle{settings.num_threads, benchmark_thread, std::ref(pq),
                                                std::ref(shared_data), instance};
    task_handle.wait();

    std::clog << "done" << std::endl;

    bool success = true;
    std::clog << "Verifying..." << std::flush;
    bool valid = true;
    if (shared_data.processed_nodes.load() + shared_data.ignored_nodes.load() != shared_data.pushed_nodes.load()) {
        std::clog << "failed: Not all nodes were popped" << std::endl;
        success = false;
        valid = false;
    } else {
        if (shared_data.best_value.load() == settings.solution) {
            std::clog << "done" << std::endl;
        } else {
            std::clog << "failed: Wrong solution" << std::endl;
            success = false;
            valid = false;
        }
    }
    std::clog << '\n';
    print_shared_data(settings, shared_data, instance, valid);
    return success;
}

void print_header() {
    std::clog << "Built on " << __DATE__ << ' ' << __TIME__ << " with:\n";
#ifdef NDEBUG
    std::clog << "  Release build\n";
#else
    std::clog << "  Debug build\n";
#endif
#ifdef __clang__
    std::clog << "  Clang " << __clang_version__ << '\n';
#elif defined(__GNUC__)
    std::clog << "  GCC " << __VERSION__ << '\n';
#else
    std::clog << "  Unknown compiler\n";
#endif
    std::clog << "  Priority queue: " << pq_name << '\n';
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
      ("j,threads", "The number of threads", cxxopts::value<int>(settings.num_threads), "NUMBER")
      ("file", "The input graph", cxxopts::value<std::filesystem::path>(settings.knapsack_file), "PATH")
      ("solution", "The reference solution", cxxopts::value<unsigned long>(settings.solution), "NUMBER")
      ("h,help", "Print this help");
    // clang-format on
    add_pq_options(options);
    options.parse_positional({"file", "solution"});

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
        pq_config = get_pq_options(result);
    }
    if (settings.knapsack_file.empty()) {
        std::cerr << "Error: No instance file specified" << std::endl;
        std::cerr << options.help() << std::endl;
        return 1;
    }

    print_settings(settings);

    bool success = run_benchmark(settings, pq_config);

    return success ? 0 : 1;
}
