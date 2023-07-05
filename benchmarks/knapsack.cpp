#include "build_info.hpp"
#include "knapsack_instance.hpp"
#include "priority_queue_selection.hpp"
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

using pq_type = PriorityQueue<unsigned long, unsigned long, false>;
using handle_type = pq_type::handle_type;

using ctx_type = thread_coordination::Context;

static_assert(sizeof(unsigned long) >= sizeof(std::uint64_t), "64bit unsigned long required");
using payload_type = unsigned long;

static constexpr std::uint8_t bits_for_index = 16;
static constexpr std::uint8_t bits_for_free_capacity = 24;
static constexpr std::uint8_t bits_for_value = 24;
static_assert(bits_for_index + bits_for_free_capacity + bits_for_value <= std::numeric_limits<payload_type>::digits,
              "Too many bits required for payload");

constexpr payload_type to_payload(std::size_t index, long long free_capacity, long long value) noexcept {
    assert(index < (1UL << bits_for_index));
    assert(free_capacity < static_cast<long long>(1UL << bits_for_free_capacity));
    assert(value < static_cast<long long>(1UL << bits_for_value));
    auto payload = static_cast<payload_type>(value);
    payload <<= bits_for_free_capacity;
    payload |= static_cast<payload_type>(free_capacity);
    payload <<= bits_for_index;
    payload |= index;
    return payload;
}

struct Settings {
    int num_threads = 4;
    std::filesystem::path knapsack_file;
    int seed = 1;
    long long solution = 0;
    pq_type::config_type pq_settings{};
};

void write_settings(Settings const& settings, std::ostream& out) {
    out << "Threads: " << settings.num_threads << '\n'
        << "Seed: " << settings.seed << '\n'
        << "Problem file: " << settings.knapsack_file << '\n'
        << "Reference solution: " << settings.solution << '\n';
}

struct ThreadStats {
    thread_coordination::time_result_type work_time{};
    long long pushed_nodes{0};
    long long ignored_nodes{0};
    long long processed_nodes{0};
};

struct SharedData {
    KnapsackInstance instance;
    std::atomic_llong solution{0};
    termination_detection::Data termination_detection_data{};

    void update_solution(long long current_solution, long long new_solution) noexcept {
        while (new_solution > current_solution &&
               !solution.compare_exchange_weak(current_solution, new_solution, std::memory_order_relaxed)) {
        }
    }
};

long long lower_bound(KnapsackInstance const& instance, long long capacity, std::size_t index) noexcept {
    long long value{0};
    while (index < instance.size() && instance.items()[index].weight <= capacity) {
        capacity -= instance.items()[index].weight;
        value += instance.items()[index].value;
        ++index;
    }
    return value;
}

long long upper_bound(KnapsackInstance const& instance, long long capacity, std::size_t index) noexcept {
    assert(index <= instance.size());
    long long value_offset = instance.prefix_sum()[index].value;
    long long target_capacity = instance.prefix_sum()[index].weight + capacity;
    while (index != instance.prefix_sum().size()) {
        if (instance.prefix_sum()[index].weight > target_capacity) {
            double fraction = static_cast<double>(target_capacity - instance.prefix_sum()[index - 1].weight) /
                static_cast<double>(instance.items()[index - 1].weight);
            return (instance.prefix_sum()[index - 1].value - value_offset) +
                static_cast<long long>(static_cast<double>(instance.items()[index - 1].value) * fraction);
        }
        ++index;
    }
    // All items fit
    return instance.prefix_sum()[index - 1].value - value_offset;
}

bool process_node(handle_type& handle, ThreadStats& stats, SharedData& data) {
    auto node = handle.try_pop();
    if (!node) {
        return false;
    }
    auto solution = data.solution.load(std::memory_order_relaxed);
    auto ub = static_cast<long long>(node->first);
    if (ub <= solution) {
        // The upper bound of this node is worse than the current solution
        ++stats.ignored_nodes;
        return true;
    }
    ++stats.processed_nodes;
    auto e = node->second;
    std::size_t index = e & ((1UL << bits_for_index) - 1);
    e >>= bits_for_index;
    auto free_capacity = static_cast<long long>(e & ((1UL << bits_for_free_capacity) - 1));
    assert(free_capacity <= data.instance.capacity());
    e >>= bits_for_free_capacity;
    auto value = static_cast<long long>(e & ((1UL << bits_for_value) - 1));
    if (index == data.instance.size()) {
        data.update_solution(solution, static_cast<long long>(value));
        return true;
    }
    if (data.instance.items()[index].weight <= free_capacity) {
        auto new_value = value + data.instance.items()[index].value;
        auto new_capacity = free_capacity - data.instance.items()[index].weight;
        auto payload = to_payload(index + 1, new_capacity, new_value);
        handle.push({ub, payload});
        ++stats.pushed_nodes;
    }
    auto upper_bound_without_next = value + upper_bound(data.instance, free_capacity, index + 1);
    if (upper_bound_without_next <= solution) {
        return true;
    }
    handle.push({upper_bound_without_next, to_payload(index + 1, free_capacity, value)});
    ++stats.pushed_nodes;
    return true;
}

ThreadStats benchmark_thread(ctx_type ctx, pq_type& pq, SharedData& data) {
    ThreadStats stats;
    auto handle = pq.get_handle();
    if (ctx.get_id() == 0) {
        auto ub = upper_bound(data.instance, data.instance.capacity(), 0);
        auto lb = lower_bound(data.instance, data.instance.capacity(), 0);
        data.solution.store(lb, std::memory_order_relaxed);
        handle.push({ub, to_payload(0, data.instance.capacity(), 0)});
        ++stats.pushed_nodes;
        stats.work_time = ctx.execute_synchronized([&]() {
            while (termination_detection::try_do(ctx.get_num_threads(), data.termination_detection_data,
                                                 [&]() { return process_node(handle, stats, data); })) {
            }
        });
    } else {
        stats.work_time = ctx.execute_synchronized([&]() {
            while (termination_detection::try_do(ctx.get_num_threads(), data.termination_detection_data,
                                                 [&]() { return process_node(handle, stats, data); })) {
            }
        });
    }
    return stats;
}

void write_stats(std::vector<ThreadStats> const& stats, std::ostream& out) {
    out << "thread,start_time,end_time,pushed,processed,ignored\n";
    for (std::size_t i = 0; i < stats.size(); ++i) {
        auto const& s = stats[i];
        // clang-format off
        out << i << ','
            << s.work_time.start.time_since_epoch().count() << ','
            << s.work_time.end.time_since_epoch().count() << ','
            << s.pushed_nodes << ','
            << s.processed_nodes << ','
            << s.ignored_nodes
            << '\n';
        // clang-format on
    }
}

bool run_benchmark(Settings const& settings) {
    SharedData shared_data;
    try {
        shared_data.instance = KnapsackInstance(settings.knapsack_file);
    } catch (std::runtime_error const& e) {
        std::clog << "Error reading problem file: " << e.what() << std::endl;
        return false;
    }

    if (shared_data.instance.size() + 1 >= (1 << bits_for_index) ||
        shared_data.instance.capacity() >= (1 << bits_for_free_capacity) ||
        settings.solution >= (1 << bits_for_value)) {
        std::clog << "Error: Instance cannot be represented\n";
        return false;
    }
    std::clog << "done\n";
    auto pq = pq_type(settings.num_threads, 0, settings.pq_settings);
    std::clog << "Priority queue: ";
    pq.describe(std::clog) << '\n';
    std::vector<ThreadStats> all_stats(static_cast<std::size_t>(settings.num_threads));
    thread_coordination::affinity::NUMA numa_affinity{CORES_PER_NUMA_NODE, NUM_NUMA_NODES};
    std::clog << "Solving..." << std::flush;
    auto t_start = std::chrono::steady_clock::now();
    thread_coordination::TaskHandle task_handle{
        numa_affinity, settings.num_threads,
        [&](auto ctx) { all_stats[static_cast<std::size_t>(ctx.get_id())] = benchmark_thread(ctx, pq, shared_data); }};
    task_handle.wait();
    auto t_end = std::chrono::steady_clock::now();
    std::clog << "done (" << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count() << "ms)"
              << std::endl;
    ThreadStats first = all_stats.front();
    auto accum_stats = std::accumulate(all_stats.begin() + 1, all_stats.end(), first, [](auto accum, auto const& e) {
        accum.work_time.start = std::min(accum.work_time.start, e.work_time.start);
        accum.work_time.end = std::max(accum.work_time.end, e.work_time.end);
        accum.pushed_nodes += e.pushed_nodes;
        accum.processed_nodes += e.processed_nodes;
        accum.ignored_nodes += e.ignored_nodes;
        return accum;
    });
    if (accum_stats.processed_nodes + accum_stats.ignored_nodes != accum_stats.pushed_nodes) {
        std::clog << "Error: Not all nodes were popped" << std::endl;
        return false;
    }
    if (shared_data.solution.load() != settings.solution) {
        std::clog << "Error: Wrong solution (got " << shared_data.solution.load() << ", expected " << settings.solution
                  << ")" << std::endl;
        return false;
    }
    write_stats(all_stats, std::cout);
    return true;
}

int main(int argc, char* argv[]) {
    write_build_info(std::clog);
    Settings settings{};
    cxxopts::Options cmd(argv[0]);
    // clang-format off
    cmd.add_options()
      ("j,threads", "The number of threads", cxxopts::value<int>(settings.num_threads), "NUMBER")
      ("file", "The input graph", cxxopts::value<std::filesystem::path>(settings.knapsack_file), "PATH")
      ("solution", "The reference solution", cxxopts::value<long long>(settings.solution), "NUMBER")
      ("h,help", "Print this help");
    // clang-format on
    pq_type::add_options(cmd, settings.pq_settings);
    cmd.parse_positional({"file", "solution"});

    try {
        auto args = cmd.parse(argc, argv);
        if (args.count("help") > 0) {
            std::cerr << cmd.help() << std::endl;
            return EXIT_SUCCESS;
        }
    } catch (cxxopts::OptionParseException const& e) {
        std::cerr << "Error parsing arguments: " << e.what() << '\n';
        std::cerr << cmd.help() << std::endl;
        return EXIT_FAILURE;
    }

    std::clog << "\nCommand line:";
    for (int i = 0; i < argc; ++i) {
        std::clog << ' ' << argv[i];
    }
    std::clog << '\n';
    write_settings(settings, std::clog);

    if (settings.knapsack_file.empty()) {
        std::cerr << "Error: No instance file specified" << std::endl;
        std::cerr << cmd.help() << std::endl;
        return EXIT_FAILURE;
    }

    bool success = run_benchmark(settings);
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
