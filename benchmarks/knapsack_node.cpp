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

static constexpr auto default_n = 1000;
static constexpr auto default_f = 0.5;
#ifdef INTEGER_INSTANCE
using data_type = long long;
static constexpr data_type default_a = 1000;
static constexpr data_type default_b = 100000;
static constexpr data_type default_l = 10000;
static constexpr data_type default_u = 12500;
#else
using data_type = double;
static constexpr data_type default_a = 0.01;
static constexpr data_type default_b = 1.01;
static constexpr data_type default_l = 0.1;
static constexpr data_type default_u = 0.125;
#endif

struct Node {
    data_type upper_bound;
    std::size_t index;
    data_type free_capacity;
    data_type value;
};

bool operator==(Node const& lhs, Node const& rhs) noexcept {
    return lhs.upper_bound == rhs.upper_bound && lhs.index == rhs.index && lhs.free_capacity == rhs.free_capacity &&
        lhs.value == rhs.value;
}

struct NodePriority {
    static auto get(Node const& node) noexcept {
        return node.upper_bound;
    }
};

using pq_type = PQWrapper<false, data_type, Node, NodePriority>;
using handle_type = pq_type::handle_type;

struct Settings {
    int num_threads = 4;
    long long n = default_n;
    data_type a = default_a;
    data_type b = default_b;
    data_type l = default_l;
    data_type u = default_u;
    double f = default_f;
    unsigned long s = 1;

    void options(cxxopts::Options& cmd) {
        cmd.add_options()("j,threads", "The number of threads", cxxopts::value<int>(num_threads), "NUMBER")(
            "n,num-elements", "Number of elements", cxxopts::value<long long>(n), "NUMBER")(
            "a,weight-min", "Min weight", cxxopts::value<data_type>(a), "NUMBER")(
            "b,weight-max", "Max weight", cxxopts::value<data_type>(b), "NUMBER")(
            "l,p-min", "Min add to profits", cxxopts::value<data_type>(l), "NUMBER")(
            "u,p-max", "Max add to profits", cxxopts::value<data_type>(u), "NUMBER")(
            "f,capacity-divide", "Fraction of weights", cxxopts::value<double>(f), "NUMBER")(
            "s,seed", "Seed", cxxopts::value<unsigned long>(s), "NUMBER");
    }
};

void write_settings(Settings const& settings, std::ostream& out) {
    out << "Threads: " << settings.num_threads << '\n'
        << "Number of items: " << settings.n << '\n'
        << "Instance parameters: "
        << "weight_min=" << settings.a << " weight_max=" << settings.b << " profit_add_min=" << settings.l
        << " profit_add_max=" << settings.u << " cap_fraction=" << settings.f << '\n'
        << "Seed: " << settings.s << '\n';
}

struct ThreadStats {
    long long pushed_nodes{0};
    long long processed_nodes{0};
    long long ignored_nodes{0};
};

struct SharedData {
    std::atomic<data_type> solution{0};
    termination_detection::Data termination_detection_data{};

    void update_solution(data_type& current, data_type update) noexcept {
        while (update > current) {
            if (solution.compare_exchange_weak(current, update, std::memory_order_relaxed)) {
                current = update;
            }
        }
    }
};

bool process_node(handle_type& handle, ThreadStats& stats, SharedData& data,
                  KnapsackInstance<data_type> const& instance) {
    auto node = handle.try_pop();
    if (!node) {
        return false;
    }
    ++stats.processed_nodes;
    auto solution = data.solution.load(std::memory_order_relaxed);
    if (node->upper_bound <= solution) {
        ++stats.ignored_nodes;
        return true;
    }
    auto [lb, ub] = instance.compute_bounds_linear(node->free_capacity, node->index + 1);
    data.update_solution(solution, node->value + lb);
    if (node->index + 2 < instance.size()) {
        if (node->value + ub > solution) {
            handle.push({node->value + ub, node->index + 1, node->free_capacity, node->value});
            ++stats.pushed_nodes;
        }
        if (node->free_capacity >= instance.weight(node->index)) {
            handle.push({node->upper_bound, node->index + 1, node->free_capacity - instance.weight(node->index),
                         node->value + instance.value(node->index)});
            ++stats.pushed_nodes;
        }
    }
    return true;
}

ThreadStats benchmark_thread(task::Control tc, pq_type& pq, SharedData& data,
                             KnapsackInstance<data_type> const& instance) {
    ThreadStats stats;
    handle_type handle = pq.get_handle();
    if (tc.id() == 0) {
        auto [lb, ub] = instance.compute_bounds_linear(instance.capacity(), 0);
        data.solution.store(lb, std::memory_order_relaxed);
        if (ub > lb) {
            handle.push({ub, 0, instance.capacity(), 0});
            ++stats.pushed_nodes;
        }
    }
    tc.synchronize();
    while (termination_detection::try_do(tc.num_threads(), data.termination_detection_data,
                                         [&]() { return process_node(handle, stats, data, instance); })) {
    }
    tc.synchronize();
    return stats;
}

bool run_benchmark(Settings const& settings, pq_type& pq, KnapsackInstance<data_type> const& instance) {
    std::vector<ThreadStats> all_stats(static_cast<std::size_t>(settings.num_threads));
    affinity::NUMA numa_affinity{cores_per_numa_node, num_numa_nodes};
    SharedData shared_data;
    std::clog << "Working...\n";
    auto start_time = std::chrono::steady_clock::now();
    task::Runner runner{numa_affinity, settings.num_threads, [&](auto tc) {
                            all_stats[static_cast<std::size_t>(tc.id())] =
                                benchmark_thread(tc, pq, shared_data, instance);
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
        std::clog << "Error: " << accum_stats.pushed_nodes - accum_stats.processed_nodes << " nodes were not processed"
                  << std::endl;
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
    cmd.add_options()("h,help", "Print this help");
    settings.options(cmd);
    add_options(cmd);

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

    auto pq = create<false, data_type, Node, NodePriority>(settings.num_threads, 1 << 24, args);
    std::clog << "Priority queue: ";
    describe(pq, std::clog) << '\n' << '\n';
    std::clog << "Generating instance...\n";
    auto instance =
        KnapsackInstance<data_type>(settings.n, settings.a, settings.b, settings.l, settings.u, settings.f, settings.s);
    bool success = run_benchmark(settings, pq, instance);
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
