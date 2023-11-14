#include "build_info.hpp"
#include "knapsack_instance.hpp"
#include "task.hpp"
#include "termination_detection.hpp"
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

using pq_type = PQWrapper<false>;
using handle_type = pq_type::handle_type;

static_assert(sizeof(unsigned long) >= sizeof(std::uint64_t), "64bit unsigned long required");
using payload_type = unsigned long;

constexpr pq_type::value_type to_payload(long long upper_bound, std::size_t index, long long free_capacity,
                                         long long value) noexcept {
    assert(upper_bound >= 0 && static_cast<unsigned long>(upper_bound) < 1UL << 32);
    assert(index < (1UL << 32));
    assert(free_capacity >= 0 && static_cast<unsigned long>(free_capacity) < 1UL << 32);
    assert(value >= 0 && static_cast<unsigned long>(value) < 1UL << 32);
    return pq_type::value_type{static_cast<unsigned long>(upper_bound) << 32 | static_cast<payload_type>(index),
                               static_cast<payload_type>(free_capacity) << 32 | static_cast<payload_type>(value)};
}

struct Settings {
    int num_threads = 4;
    std::filesystem::path knapsack_file;
    int seed = 1;
};

void write_settings(Settings const& settings, std::ostream& out) {
    out << "Threads: " << settings.num_threads << '\n'
        << "Seed: " << settings.seed << '\n'
        << "Problem file: " << settings.knapsack_file << '\n';
}

struct ThreadStats {
    long long pushed_nodes{0};
    long long processed_nodes{0};
    long long ignored_nodes{0};
};

struct SharedData {
    KnapsackInstance<long long> instance;
    std::atomic_llong solution{0};
    termination_detection::Data termination_detection_data{};

    void update_solution(long long& current, long long update) noexcept {
        while (update > current) {
            if (solution.compare_exchange_weak(current, update, std::memory_order_relaxed)) {
                current = update;
            }
        }
    }
};

bool process_node(handle_type& handle, ThreadStats& stats, SharedData& data) {
    auto node = handle.try_pop();
    if (!node) {
        return false;
    }
    ++stats.processed_nodes;
    auto solution = data.solution.load(std::memory_order_relaxed);
    auto upper_bound = static_cast<long long>(node->first >> 32);
    if (upper_bound <= solution) {
        ++stats.ignored_nodes;
        return true;
    }
    auto index = static_cast<std::size_t>(node->first & ((1UL << 32) - 1));
    auto free_capacity = static_cast<long long>(node->second >> 32);
    assert(free_capacity <= data.instance.capacity());
    auto value = static_cast<long long>(node->second & ((1UL << 32) - 1));
    auto [lb, ub] = data.instance.compute_bounds_linear(free_capacity, index + 1);
    data.update_solution(solution, value + lb);
    if (index + 2 < data.instance.size()) {
        if (value + ub > solution) {
            handle.push(to_payload(value + ub, index + 1, free_capacity, value));
            ++stats.pushed_nodes;
        }
        if (free_capacity >= data.instance.weight(index)) {
            handle.push(to_payload(upper_bound, index + 1, free_capacity - data.instance.weight(index),
                                   value + data.instance.value(index)));
            ++stats.pushed_nodes;
        }
    }
    return true;
}

ThreadStats benchmark_thread(task::Control tc, pq_type& pq, SharedData& data) {
    ThreadStats stats;
    handle_type handle = pq.get_handle();
    if (tc.id() == 0) {
        auto [lb, ub] = data.instance.compute_bounds_linear(data.instance.capacity(), 0);
        data.solution.store(lb, std::memory_order_relaxed);
        if (ub > lb) {
            handle.push(to_payload(ub, 0, data.instance.capacity(), 0));
            ++stats.pushed_nodes;
        }
    }
    tc.synchronize();
    while (termination_detection::try_do(tc.num_threads(), data.termination_detection_data,
                                         [&]() { return process_node(handle, stats, data); })) {
    }
    tc.synchronize();
    return stats;
}

bool run_benchmark(Settings const& settings, pq_type& pq) {
    SharedData shared_data;
    try {
        shared_data.instance = KnapsackInstance(settings.knapsack_file);
    } catch (std::runtime_error const& e) {
        std::clog << "Error reading problem file: " << e.what() << std::endl;
        return false;
    }

    if (shared_data.instance.size() + 1 >= (1UL << 32) || shared_data.instance.capacity() >= (1LL << 32)) {
        std::clog << "Error: Instance cannot be represented\n";
        return false;
    }
    std::vector<ThreadStats> all_stats(static_cast<std::size_t>(settings.num_threads));
    affinity::NUMA numa_affinity{cores_per_numa_node, num_numa_nodes};
    std::clog << "Working...\n";
    auto start_time = std::chrono::steady_clock::now();
    task::Runner runner{numa_affinity, settings.num_threads, [&](auto tc) {
                            all_stats[static_cast<std::size_t>(tc.id())] = benchmark_thread(tc, pq, shared_data);
                        }};
    runner.wait();
    auto end_time = std::chrono::steady_clock::now();
    std::clog << "Finished\n" << std::endl;
    auto accum_stats =
        std::accumulate(all_stats.begin() + 1, all_stats.end(), all_stats.front(), [](auto accum, auto const& e) {
            accum.pushed_nodes += e.pushed_nodes;
            accum.processed_nodes += e.processed_nodes;
            accum.ignored_nodes += e.ignored_nodes;
            return accum;
        });
    std::clog << "Time (s): " << std::fixed << std::setprecision(3)
              << std::chrono::duration<double>(end_time - start_time).count() << '\n';
    std::clog << "Solution: " << shared_data.solution.load() << '\n';
    std::clog << "Processed nodes: " << accum_stats.processed_nodes << '\n';
    std::clog << "Ignored nodes: " << accum_stats.ignored_nodes << '\n';
    if (accum_stats.processed_nodes != accum_stats.pushed_nodes) {
        std::clog << "Error: Not all nodes were popped" << std::endl;
        return false;
    }
    std::cout << "time,processed,ignored,solution\n";
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << ','
              << accum_stats.processed_nodes << ',' << accum_stats.ignored_nodes << ',' << shared_data.solution.load()
              << '\n';
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
    add_options(cmd);
    cmd.parse_positional({"file"});

    auto args = cxxopts::ParseResult{};
    try {
        args = cmd.parse(argc, argv);
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

    auto pq = create<false>(settings.num_threads, 1 << 24, args);
    std::clog << "Priority queue: ";
    describe(pq, std::clog) << '\n' << '\n';
    bool success = run_benchmark(settings, pq);
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
