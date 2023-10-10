#include "build_info.hpp"
#include "knapsack_instance.hpp"
#include "task.hpp"
#include "termination_detection.hpp"
#include "wrapper/priority.hpp"
#include "wrapper/selector.hpp"

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

#ifdef CORES_PER_NUMA_NODE
static constexpr auto cores_per_numa_node = CORES_PER_NUMA_NODE;
#else
static constexpr auto cores_per_numa_node = 4;
#endif
#ifdef NUM_NUMA_NODES
static constexpr auto num_numa_nodes = NUM_NUMA_NODES;
#else
static constexpr auto num_numa_nodes = 16;
#endif
using pq_type = PriorityQueue<unsigned long, unsigned long, Priority::Max>;
using handle_type = pq_type::handle_type;

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
    pq_type::config_type pq_settings{};
};

void write_settings(Settings const& settings, std::ostream& out) {
    out << "Threads: " << settings.num_threads << '\n'
        << "Seed: " << settings.seed << '\n'
        << "Problem file: " << settings.knapsack_file << '\n';
}

struct ThreadStats {
    std::pair<std::chrono::high_resolution_clock::time_point, std::chrono::high_resolution_clock::time_point> work_time;
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
    auto upper_bound_without_next = value + data.instance.upper_bound(free_capacity, index + 1);
    if (upper_bound_without_next <= solution) {
        return true;
    }
    handle.push({upper_bound_without_next, to_payload(index + 1, free_capacity, value)});
    ++stats.pushed_nodes;
    return true;
}

ThreadStats benchmark_thread(task::Control tc, pq_type& pq, SharedData& data) {
    ThreadStats stats;
    auto handle = pq.get_handle();
    if (tc.id() == 0) {
        auto ub = data.instance.upper_bound(data.instance.capacity(), 0);
        auto lb = data.instance.lower_bound(data.instance.capacity(), 0);
        data.solution.store(lb, std::memory_order_relaxed);
        handle.push({ub, to_payload(0, data.instance.capacity(), 0)});
        ++stats.pushed_nodes;
    }
    tc.synchronize();
    auto start_time = std::chrono::high_resolution_clock::now();
    while (termination_detection::try_do(tc.num_threads(), data.termination_detection_data,
                                         [&]() { return process_node(handle, stats, data); })) {
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    tc.synchronize();
    stats.work_time = {start_time, end_time};
    return stats;
}

void write_stats(ThreadStats const& stats, long long solution, std::ostream& out) {
    out << "solution,time,pushed,processed,ignored\n";
    out << solution << ',' << (stats.work_time.second - stats.work_time.first).count() << ',' << stats.pushed_nodes
        << ',' << stats.processed_nodes << ',' << stats.ignored_nodes << '\n';
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
        shared_data.instance.capacity() >= (1 << bits_for_free_capacity)) {
        std::clog << "Error: Instance cannot be represented\n";
        return false;
    }
    auto pq = pq_type(settings.num_threads, 1 << 20, settings.pq_settings);
    std::clog << "Priority queue: ";
    pq.describe(std::clog) << '\n' << '\n';
    std::vector<ThreadStats> all_stats(static_cast<std::size_t>(settings.num_threads));
    affinity::NUMA numa_affinity{cores_per_numa_node, num_numa_nodes};
    std::clog << "Solving...\n";
    task::Runner runner{numa_affinity, settings.num_threads, [&](auto tc) {
                            all_stats[static_cast<std::size_t>(tc.id())] = benchmark_thread(tc, pq, shared_data);
                        }};
    runner.wait();
    std::clog << "Finished\n" << std::endl;
    auto accum_stats =
        std::accumulate(all_stats.begin() + 1, all_stats.end(), all_stats.front(), [](auto accum, auto const& e) {
            accum.work_time.first = std::min(accum.work_time.first, e.work_time.first);
            accum.work_time.second = std::max(accum.work_time.second, e.work_time.second);
            accum.pushed_nodes += e.pushed_nodes;
            accum.processed_nodes += e.processed_nodes;
            accum.ignored_nodes += e.ignored_nodes;
            return accum;
        });
    std::clog << "Time (s): " << std::setprecision(3)
              << std::chrono::duration<double>(accum_stats.work_time.second - accum_stats.work_time.first).count()
              << '\n';
    std::clog << "Pushed nodes: " << accum_stats.pushed_nodes << '\n';
    std::clog << "Processed nodes: " << accum_stats.processed_nodes << '\n';
    std::clog << "Ignored nodes: " << accum_stats.ignored_nodes << '\n';
    if (accum_stats.processed_nodes + accum_stats.ignored_nodes != accum_stats.pushed_nodes) {
        std::clog << "Error: Not all nodes were popped" << std::endl;
        return false;
    }
    write_stats(accum_stats, shared_data.solution.load(), std::cout);
    return true;
}

int main(int argc, char* argv[]) {
    write_build_info(std::clog);
    std::clog << "\nCommand line:";
    for (int i = 0; i < argc; ++i) {
        std::clog << ' ' << argv[i];
    }
    std::clog << '\n' << '\n';

    Settings settings{};
    cxxopts::Options cmd(argv[0]);
    // clang-format off
    cmd.add_options()
      ("j,threads", "The number of threads", cxxopts::value<int>(settings.num_threads), "NUMBER")
      ("file", "The input graph", cxxopts::value<std::filesystem::path>(settings.knapsack_file), "PATH")
      ("h,help", "Print this help");
    // clang-format on
    pq_type::add_options(cmd, settings.pq_settings);
    cmd.parse_positional({"file"});

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

    write_settings(settings, std::clog);

    if (settings.knapsack_file.empty()) {
        std::cerr << "Error: No instance file specified" << std::endl;
        std::cerr << cmd.help() << std::endl;
        return EXIT_FAILURE;
    }

    bool success = run_benchmark(settings);
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
