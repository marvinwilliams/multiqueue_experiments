#include "util/build_info.hpp"
#include "util/knapsack_instance.hpp"
#include "util/termination_detection.hpp"
#include "util/thread_coordination.hpp"
#include "wrapper/selector.hpp"

#include <cxxopts.hpp>

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

#ifdef FLOAT_INSTANCE
using data_type = double;
static constexpr data_type default_min_weight = 0.01;
static constexpr data_type default_max_weight = 1.01;
static constexpr data_type default_min_add = 0.1;
static constexpr data_type default_max_add = 0.125;
struct Payload {
    std::size_t index;
    data_type free_capacity;
    data_type weight;
};
using pq_type = PQ<false, double, Payload>;
using node_type = pq_type::value_type;

constexpr auto to_payload(data_type upper_bound, std::size_t index, data_type free_capacity, data_type value) noexcept {
    return node_type{upper_bound, Payload{index, free_capacity, value}};
}

data_type extract_upper_bound(node_type const& node) noexcept {
    return node.first;
}

std::size_t extract_index(node_type const& node) noexcept {
    return node.second.index;
}

data_type extract_free_capacity(node_type const& node) noexcept {
    return node.second.free_capacity;
}

data_type extract_value(node_type const& node) noexcept {
    return node.second.weight;
}

#else
using data_type = unsigned long;
static constexpr data_type default_min_weight = 1000;
static constexpr data_type default_max_weight = 100000;
static constexpr data_type default_min_add = 10000;
static constexpr data_type default_max_add = 12500;
using pq_type = PQ<false, unsigned long, unsigned long>;
using node_type = pq_type::value_type;

constexpr auto to_payload(data_type upper_bound, std::size_t index, data_type free_capacity, data_type value) noexcept {
    static_assert(sizeof(data_type) >= sizeof(std::uint64_t), "64bit data_type required");
    assert(upper_bound <= std::numeric_limits<std::uint32_t>::max());
    assert(index <= std::numeric_limits<std::uint32_t>::max());
    assert(free_capacity <= std::numeric_limits<std::uint32_t>::max());
    assert(value <= std::numeric_limits<std::uint32_t>::max());

    return node_type{static_cast<std::uint64_t>(index) | (static_cast<std::uint64_t>(upper_bound) << 32),
                     static_cast<std::uint64_t>(value) | (static_cast<std::uint64_t>(free_capacity) << 32)};
}

data_type extract_upper_bound(node_type const& node) noexcept {
    return node.first >> 32;
}

std::size_t extract_index(node_type const& node) noexcept {
    return node.first & ((1UL << 32) - 1);
}

data_type extract_free_capacity(node_type const& node) noexcept {
    return node.second >> 32;
}

data_type extract_value(node_type const& node) noexcept {
    return node.second & ((1UL << 32) - 1);
}

#endif

using handle_type = pq_type::handle_type;

struct Settings {
    int num_threads = 4;
    long long n = 1000;
    data_type min_weight = default_min_weight;
    data_type max_weight = default_max_weight;
    data_type min_add = default_min_add;
    data_type max_add = default_max_add;
    double capacity_factor = 0.5;
    unsigned int seed = 1;
    pq_type::settings_type pq_settings{};
};

void register_cmd_options(Settings& settings, cxxopts::Options& cmd) {
    cmd.add_options()
        // clang-format off
        ("j,threads", "The number of threads", cxxopts::value<int>(settings.num_threads), "NUMBER")
        ("n,num-elements", "Number of elements", cxxopts::value<long long>(settings.n), "NUMBER")
        ("a,min-weight", "Min weight", cxxopts::value<data_type>(settings.min_weight), "NUMBER")
        ("b,max-weight", "Max weight", cxxopts::value<data_type>(settings.max_weight), "NUMBER")
        ("l,min-add", "Min add to profits", cxxopts::value<data_type>(settings.min_add), "NUMBER")
        ("u,max-add", "Max add to profits", cxxopts::value<data_type>(settings.max_add), "NUMBER")
        ("f,factor", "Capacity as factor of expected total weight", cxxopts::value<double>(settings.capacity_factor), "NUMBER")
        ("s,seed", "Seed", cxxopts::value<unsigned int>(settings.seed), "NUMBER");
    // clang-format on
    settings.pq_settings.register_cmd_options(cmd);
}

void write_settings_human_readable(Settings const& settings, std::ostream& out) {
    out << "Threads: " << settings.num_threads << '\n';
    out << "Number of items: " << settings.n << '\n';
    out << "Weights: [" << settings.min_weight << ", " << settings.max_weight << "]\n";
    out << "Add to profit: [" << settings.min_add << ", " << settings.max_add << "]\n";
    out << "Capacity factor: " << settings.capacity_factor << '\n';
    out << "Seed: " << settings.seed << '\n';
    settings.pq_settings.write_human_readable(out);
}

void write_settings_json(Settings const& settings, std::ostream& out) {
    out << '{';
    out << std::quoted("num-threads") << ':' << settings.num_threads << ',';
    out << std::quoted("num-elements") << ':' << settings.n << ',';
    out << std::quoted("min-weight") << ':' << settings.min_weight << ',';
    out << std::quoted("max-weight") << ':' << settings.max_weight << ',';
    out << std::quoted("min-add") << ':' << settings.min_add << ',';
    out << std::quoted("max-add") << ':' << settings.max_add << ',';
    out << std::quoted("capacity-factor") << ':' << settings.capacity_factor << ',';
    out << std::quoted("seed") << ':' << settings.seed << ',';
    out << std::quoted("pq") << ':';
    settings.pq_settings.write_json(out);
    out << '}';
}

struct Counter {
    long long pushed_nodes{0};
    long long processed_nodes{0};
    long long ignored_nodes{0};
};

struct SharedData {
    std::atomic<data_type> solution{0};
    termination_detection::TerminationDetection termination_detection;
};

void process_node(node_type const& node, handle_type& handle, Counter& counter, SharedData& data,
                  KnapsackInstance<data_type> const& instance) {
    ++counter.processed_nodes;
    auto solution = data.solution.load(std::memory_order_relaxed);
    auto upper_bound = extract_upper_bound(node);
    if (upper_bound <= solution) {
        ++counter.ignored_nodes;
        return;
    }
    auto index = extract_index(node);
    auto free_capacity = extract_free_capacity(node);
    assert(free_capacity <= instance.capacity());
    auto value = extract_value(node);
    auto [lb, ub] = instance.compute_bounds_linear(free_capacity, index + 1);
    while (value + lb > solution &&
           !data.solution.compare_exchange_weak(solution, value + lb, std::memory_order_relaxed)) {
    }
    solution = std::max(solution, value + lb);
    if (index + 2 < instance.size()) {
        if (value + ub > solution) {
            handle.push(to_payload(value + ub, index + 1, free_capacity, value));
            ++counter.pushed_nodes;
        }
        if (free_capacity >= instance.weight(index)) {
            handle.push(to_payload(upper_bound, index + 1, free_capacity - instance.weight(index),
                                   value + instance.value(index)));
            ++counter.pushed_nodes;
        }
    }
}

Counter benchmark_thread(thread_coordination::Context& thread_context, pq_type& pq, SharedData& data,
                         KnapsackInstance<data_type> const& instance) {
    Counter counter;
    handle_type handle = pq.get_handle();
    if (thread_context.id() == 0) {
        auto [lb, ub] = instance.compute_bounds_linear(instance.capacity(), 0);
        data.solution.store(lb, std::memory_order_relaxed);
        handle.push(to_payload(ub, 0, instance.capacity(), 0));
        ++counter.pushed_nodes;
    }
    std::optional<node_type> node;
    thread_context.synchronize();
    while (data.termination_detection.repeat([&]() {
        node = handle.try_pop();
        return node.has_value();
    })) {
        process_node(*node, handle, counter, data, instance);
    }
    thread_context.synchronize();
    return counter;
}

void run_benchmark(Settings const& settings) {
    std::clog << "Generating instance...\n";
    auto instance = KnapsackInstance<data_type>(settings.n, settings.min_weight, settings.max_weight, settings.min_add,
                                                settings.max_add, settings.capacity_factor, settings.seed);
    std::vector<Counter> thread_counter(static_cast<std::size_t>(settings.num_threads));
    SharedData shared_data{0, termination_detection::TerminationDetection(settings.num_threads)};
    auto pq = pq_type(settings.num_threads, std::size_t(10'000'000), settings.pq_settings);
    std::clog << "Working...\n";
    auto start_time = std::chrono::steady_clock::now();
    thread_coordination::Dispatcher dispatcher{settings.num_threads, [&](auto ctx) {
                                                   auto t_id = static_cast<std::size_t>(ctx.id());
                                                   thread_counter[t_id] =
                                                       benchmark_thread(ctx, pq, shared_data, instance);
                                               }};
    dispatcher.wait();
    auto end_time = std::chrono::steady_clock::now();
    std::clog << "Done\n";
    auto total_counts =
        std::accumulate(thread_counter.begin(), thread_counter.end(), Counter{}, [](auto sum, auto const& counter) {
            sum.pushed_nodes += counter.pushed_nodes;
            sum.processed_nodes += counter.processed_nodes;
            sum.ignored_nodes += counter.ignored_nodes;
            return sum;
        });
    std::clog << '\n';
    std::clog << "Time (s): " << std::fixed << std::setprecision(3)
              << std::chrono::duration<double>(end_time - start_time).count() << '\n';
    std::clog << "Solution: " << shared_data.solution.load() << '\n';
    std::clog << "Processed nodes: " << total_counts.processed_nodes << '\n';
    std::clog << "Ignored nodes: " << total_counts.ignored_nodes << '\n';
    if (total_counts.processed_nodes != total_counts.pushed_nodes) {
        std::cerr << "Error: Not all nodes were processed (" << total_counts.pushed_nodes - total_counts.processed_nodes
                  << ')' << '\n';
        std::exit(EXIT_FAILURE);
    }
    std::cout << '{';
    std::cout << std::quoted("settings") << ':';
    write_settings_json(settings, std::cout);
    std::cout << ',';
    std::cout << std::quoted("time-ns") << ':' << std::chrono::nanoseconds{end_time - start_time}.count() << ',';
    std::cout << std::quoted("processed-nodes") << ':' << total_counts.processed_nodes << ',';
    std::cout << std::quoted("ignored-nodes") << ':' << total_counts.ignored_nodes << ',';
    std::cout << std::quoted("solution") << ':' << shared_data.solution.load();
    std::cout << '}' << '\n';
}

int main(int argc, char* argv[]) {
    write_build_info(std::clog);
    std::clog << '\n';

    std::clog << "= Priority queue =\n";
    pq_type::write_human_readable(std::clog);
    std::clog << '\n';

    std::clog << "= Command line =\n";
    for (int i = 0; i < argc; ++i) {
        std::clog << argv[i];
        if (i != argc - 1) {
            std::clog << ' ';
        }
    }
    std::clog << '\n' << '\n';

    cxxopts::Options cmd(argv[0]);
    cmd.add_options()("h,help", "Print this help");
    Settings settings{};
    register_cmd_options(settings, cmd);

    try {
        auto args = cmd.parse(argc, argv);
        if (args.count("help") > 0) {
            std::cerr << cmd.help() << std::endl;
            return EXIT_SUCCESS;
        }
    } catch (cxxopts::OptionParseException const& e) {
        std::cerr << "Error parsing command line: " << e.what() << '\n';
        std::cerr << "Use --help for usage information" << '\n';
        return EXIT_FAILURE;
    }

    std::clog << "= Settings =\n";
    write_settings_human_readable(settings, std::clog);
    std::clog << '\n';

    std::clog << "= Running benchmark =\n";
    run_benchmark(settings);
    return EXIT_SUCCESS;
}
