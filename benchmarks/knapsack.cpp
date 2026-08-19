#include "util/benchmark.hpp"
#include "util/knapsack_instance.hpp"
#include "util/selector.hpp"
#include "util/termination_detection.hpp"
#include "util/thread_coordination.hpp"

#include <cxxopts.hpp>

#include <atomic>
#include <cassert>
#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
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

struct Settings : benchmark::CommonSettings {
    std::filesystem::path instance_file{};
    pq_type::settings_type pq_settings{};
};

void register_cmd_options(Settings& settings, cxxopts::Options& cmd) {
    settings.CommonSettings::register_cmd_options(cmd);
    // clang-format off
    cmd.add_options()
        ("instance", "Instance file", cxxopts::value<std::filesystem::path>(settings.instance_file), "FILE");
    // clang-format on
    settings.pq_settings.register_cmd_options(cmd);
    cmd.parse_positional({"instance"});
}

bool validate_settings(Settings const& settings) {
    if (!settings.CommonSettings::validate()) {
        return false;
    }
    if (settings.instance_file.empty()) {
        std::cerr << "Error: No instance file specified\n";
        return false;
    }
    return settings.pq_settings.validate();
}

void write_settings_human_readable(Settings const& settings, std::ostream& out) {
    settings.CommonSettings::write_human_readable(out);
    out << "Instance file: " << settings.instance_file << '\n';
    settings.pq_settings.write_human_readable(out);
}

void write_settings_json(Settings const& settings, json::Object& obj) {
    settings.CommonSettings::write_json(obj);
    obj.entry("instance_file", settings.instance_file);
    obj.raw("pq", [&settings](std::ostream& out) { settings.pq_settings.write_json(out); });
}

struct Counter {
    long long pushed_nodes{0};
    long long processed_nodes{0};
    long long ignored_nodes{0};
};

struct SharedData {
    KnapsackInstance<data_type> instance;
    std::atomic<data_type> solution{0};
    termination_detection::TerminationDetection termination_detection;
    std::atomic_llong missing_nodes{0};
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
            if (handle.push(to_payload(value + ub, index + 1, free_capacity, value))) {
                ++counter.pushed_nodes;
            }
        }
        if (free_capacity >= data.instance.weight(index)) {
            if (handle.push(to_payload(upper_bound, index + 1, free_capacity - data.instance.weight(index),
                                       value + data.instance.value(index)))) {
                ++counter.pushed_nodes;
            }
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
    thread_context.synchronize();
    while (true) {
        std::optional<node_type> node;
        while (data.termination_detection.repeat([&]() {
            node = handle.try_pop();
            return node.has_value();
        })) {
            process_node(*node, handle, counter, data);
        }
        data.missing_nodes.fetch_add(counter.pushed_nodes - counter.processed_nodes - counter.ignored_nodes,
                                     std::memory_order_relaxed);
        thread_context.synchronize();
        if (data.missing_nodes.load(std::memory_order_relaxed) == 0) {
            break;
        }
        thread_context.synchronize();
        if (thread_context.id() == 0) {
            data.missing_nodes.store(0, std::memory_order_relaxed);
            data.termination_detection.reset();
        }
        thread_context.synchronize();
    }
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
    thread_coordination::dispatch(settings.affinity, settings.num_threads, [&](auto ctx) {
        auto t_id = static_cast<std::size_t>(ctx.id());
        thread_counter[t_id] = benchmark_thread(ctx, pq, shared_data);
    });
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
    if (total_counts.processed_nodes + total_counts.ignored_nodes != total_counts.pushed_nodes) {
        std::cerr << "Warning: Not all nodes were popped\n";
        std::cerr << "Probably the priority queue discards duplicate keys\n";
    }
    // std::fixed matters for the floating point instance type; the integral
    // fields written before it are unaffected
    std::cout << std::fixed;
    {
        json::Object root{std::cout};
        root.object("settings", [&settings](json::Object& obj) { write_settings_json(settings, obj); });
        root.object("instance", [&shared_data](json::Object& instance) {
            instance.entry("num_items", shared_data.instance.size());
            instance.entry("capacity", shared_data.instance.capacity());
        });
        root.object("results", [&](json::Object& results) {
            results.entry("time_ns", std::chrono::nanoseconds{end_time - start_time}.count());
            results.entry("processed_nodes", total_counts.processed_nodes);
            results.entry("ignored_nodes", total_counts.ignored_nodes);
            results.entry("solution", shared_data.solution.load());
        });
    }
    std::cout << '\n';
}

int main(int argc, char* argv[]) {
    benchmark::write_run_header<pq_type>(argc, argv, std::clog);

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

    if (!validate_settings(settings)) {
        return EXIT_FAILURE;
    }

    std::clog << "= Running benchmark =\n";
    run_benchmark(settings);
    return EXIT_SUCCESS;
}
