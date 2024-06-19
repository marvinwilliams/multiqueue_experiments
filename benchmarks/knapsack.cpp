#include "util/build_info.hpp"
#include "util/knapsack_instance.hpp"
#include "util/selector.hpp"
#include "util/termination_detection.hpp"
#include "util/thread_coordination.hpp"

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
struct Payload {
    std::size_t index;
    data_type free_capacity;
    data_type weight;
    friend bool operator==(Payload const& lhs, Payload const& rhs) noexcept {
        return lhs.index == rhs.index && lhs.free_capacity == rhs.free_capacity && lhs.weight == rhs.weight;
    }
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
    std::filesystem::path instance_file{};
    pq_type::settings_type pq_settings{};
};

void register_cmd_options(Settings& settings, cxxopts::Options& cmd) {
    // clang-format off
    cmd.add_options()
        ("j,threads", "The number of threads", cxxopts::value<int>(settings.num_threads), "NUMBER")
        ("instance", "Instance file", cxxopts::value<std::filesystem::path>(settings.instance_file), "FILE");
    // clang-format on
    settings.pq_settings.register_cmd_options(cmd);
    cmd.parse_positional({"instance"});
}

void write_settings_human_readable(Settings const& settings, std::ostream& out) {
    out << "Threads: " << settings.num_threads << '\n';
    out << "Instance file: " << settings.instance_file << '\n';
    settings.pq_settings.write_human_readable(out);
}

void write_settings_json(Settings const& settings, std::ostream& out) {
    out << '{';
    out << std::quoted("instance_file") << ':' << settings.instance_file << ',';
    out << std::quoted("pq") << ':';
    settings.pq_settings.write_json(out);
    out << '}';
}

struct Counter {
    long long pushed_nodes{0};
    long long processed_nodes{0};
    long long ignored_nodes{0};
    std::chrono::nanoseconds total_time{0};
    std::chrono::nanoseconds pq_time{0};
};

struct SharedData {
    KnapsackInstance<data_type> instance;
    std::atomic<data_type> solution{0};
    termination_detection::TerminationDetection termination_detection;
};

void process_node(node_type const& node, handle_type& handle, Counter& counter, SharedData& data) {
    auto solution = data.solution.load(std::memory_order_relaxed);
    auto upper_bound = extract_upper_bound(node);
    if (upper_bound <= solution) {
        ++counter.ignored_nodes;
        return;
    }
    auto index = extract_index(node);
    auto free_capacity = extract_free_capacity(node);
    assert(free_capacity <= data.instance.capacity());
    auto value = extract_value(node);
    auto [lb, ub] = data.instance.compute_bounds_linear(free_capacity, index + 1);
    while (value + lb > solution) {
        if (data.solution.compare_exchange_weak(solution, value + lb, std::memory_order_relaxed)) {
            solution = value + lb;
            break;
        }
    }
    if (index + 2 < data.instance.size()) {
        if (value + ub > solution) {
            auto now = std::chrono::steady_clock::now();
            handle.push(to_payload(value + ub, index + 1, free_capacity, value));
            counter.pq_time += std::chrono::steady_clock::now() - now;
            ++counter.pushed_nodes;
        }
        if (free_capacity >= data.instance.weight(index)) {
            auto now = std::chrono::steady_clock::now();
            handle.push(to_payload(upper_bound, index + 1, free_capacity - data.instance.weight(index),
                                   value + data.instance.value(index)));
            counter.pq_time += std::chrono::steady_clock::now() - now;
            ++counter.pushed_nodes;
        }
    }
    ++counter.processed_nodes;
}

[[gnu::noinline]] Counter benchmark_thread(thread_coordination::Context& thread_context, pq_type& pq,
                                           SharedData& data) {
    Counter counter;
    handle_type handle = pq.get_handle();
    if (thread_context.id() == 0) {
        auto [lb, ub] = data.instance.compute_bounds_linear(data.instance.capacity(), 0);
        data.solution.store(lb, std::memory_order_relaxed);
        handle.push(to_payload(ub, 0, data.instance.capacity(), 0));
        ++counter.pushed_nodes;
    }
    std::optional<node_type> node;
    thread_context.synchronize();
    while (data.termination_detection.repeat([&]() {
        auto start = std::chrono::steady_clock::now();
        node = handle.try_pop();
        auto now = std::chrono::steady_clock::now();
        if (node.has_value()) {
            counter.pq_time += now - start;
        }
        return node.has_value();
    })) {
        process_node(*node, handle, counter, data);
    }
    thread_context.synchronize();
    return counter;
}

void run_benchmark(Settings const& settings) {
    KnapsackInstance<data_type> instance;
    std::clog << "Reading instance...\n";
    try {
        instance = KnapsackInstance<data_type>(settings.instance_file);
    } catch (std::exception const& e) {
        std::cerr << "Error reading instance file: " << e.what() << '\n';
        std::exit(EXIT_FAILURE);
    }
    std::clog << "Instance has " << instance.size() << " items and " << std::fixed << instance.capacity()
              << " capacity\n";
    SharedData shared_data{std::move(instance), 0, termination_detection::TerminationDetection(settings.num_threads)};
    std::vector<Counter> thread_counter(static_cast<std::size_t>(settings.num_threads));
    auto pq = pq_type(settings.num_threads, std::size_t(10'000'000), settings.pq_settings);
    std::clog << "Working...\n";
    auto start_time = std::chrono::steady_clock::now();
    thread_coordination::Dispatcher dispatcher{settings.num_threads, [&](auto ctx) {
                                                   auto t_id = static_cast<std::size_t>(ctx.id());
                                                   thread_counter[t_id] = benchmark_thread(ctx, pq, shared_data);
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
    std::clog << "= Results =\n";
    std::clog << "Time (s): " << std::fixed << std::setprecision(3)
              << std::chrono::duration<double>(end_time - start_time).count() << '\n';
    std::clog << "Solution: " << shared_data.solution.load() << '\n';
    std::clog << "Processed nodes: " << total_counts.processed_nodes << '\n';
    std::clog << "Ignored nodes: " << total_counts.ignored_nodes << '\n';
    std::clog << "Time in priority queue:";
    for (auto const& counter : thread_counter) {
        std::clog << ' ' << std::chrono::duration<double>(counter.pq_time).count();
    }
    std::clog << '\n';
    std::clog << "Total time in priority queue: " << std::chrono::duration<double>(total_counts.pq_time).count()
              << '\n';
    std::clog << "Total time per thread:";
    for (auto const& counter : thread_counter) {
        std::clog << ' ' << std::chrono::duration<double>(counter.total_time).count();
    }
    std::clog << '\n';
    std::clog << "Total time: " << std::chrono::duration<double>(total_counts.total_time).count() << '\n';
    std::clog << "Percentage in priority queue:";
    for (auto const& counter : thread_counter) {
        std::clog << ' ' << std::fixed << std::setprecision(2)
                  << 100.0 * counter.pq_time.count() / counter.total_time.count();
    }
    std::clog << '\n';
    std::clog << "Average percentage in priority queue: " << std::fixed << std::setprecision(2)
              << 100.0 * total_counts.pq_time.count() / total_counts.total_time.count() << '\n';
    if (total_counts.processed_nodes + total_counts.ignored_nodes != total_counts.pushed_nodes) {
        std::cerr << "Warning: Not all nodes were popped\n";
        std::cerr << "Probably the priority queue discards duplicate keys\n";
    }
    std::cout << '{';
    std::cout << std::quoted("settings") << ':';
    write_settings_json(settings, std::cout);
    std::cout << ',';
    std::cout << std::quoted("instance") << ':';
    std::cout << '{';
    std::cout << std::quoted("num_items") << ':' << shared_data.instance.size() << ',';
    std::cout << std::quoted("capacity") << ':' << std::fixed << shared_data.instance.capacity();
    std::cout << '}' << ',';
    std::cout << std::quoted("results") << ':';
    std::cout << '{';
    std::cout << std::quoted("time_ns") << ':' << std::chrono::nanoseconds{end_time - start_time}.count() << ',';
    std::cout << std::quoted("processed_nodes") << ':' << total_counts.processed_nodes << ',';
    std::cout << std::quoted("ignored_nodes") << ':' << total_counts.ignored_nodes << ',';
    std::cout << std::quoted("solution") << ':' << shared_data.solution.load();
    std::cout << '}';
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
            std::cerr << cmd.help() << '\n';
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
    if (settings.instance_file.empty()) {
        std::cerr << "Error: No instance file specified" << '\n';
        std::cerr << "Use --help for usage information" << '\n';
        return EXIT_FAILURE;
    }
    std::clog << "= Running benchmark =\n";
    run_benchmark(settings);
    return EXIT_SUCCESS;
}
