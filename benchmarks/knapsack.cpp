#include "util/benchmark.hpp"
#include "util/knapsack_instance.hpp"
#include "util/memory_stats.hpp"
#include "util/parallel_search.hpp"
#include "util/thread_coordination.hpp"
#include "wrapper/selector.hpp"

#include <cxxopts.hpp>

#include <atomic>
#include <cassert>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

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

using handle_type = pq_type::handle_type;

struct Settings {
    benchmark::BaseSettings base_settings{};
    pq_type::settings_type pq_settings{};
    std::filesystem::path instance_file{};

    void register_cmd_options(cxxopts::Options& cmd) {
        base_settings.register_cmd_options(cmd);
        pq_settings.register_cmd_options(cmd);
        // clang-format off
        cmd.add_options()
            ("instance", "Instance file", cxxopts::value<std::filesystem::path>(instance_file), "FILE");
        // clang-format on
        cmd.parse_positional({"instance"});
    }

    bool validate() const {
        if (!base_settings.validate()) {
            return false;
        }
        if (!pq_settings.validate()) {
            return false;
        }
        if (instance_file.empty()) {
            std::cerr << "Error: No instance file specified\n";
            return false;
        }
        return true;
    }

    void write_human_readable(std::ostream& out) const {
        base_settings.write_human_readable(out);
        pq_settings.write_human_readable(out);
        out << "Instance file: " << instance_file << '\n';
    }

    void write_json(json::Object& obj) const {
        base_settings.write_json(obj);
        obj.object("pq", [this](json::Object& o) { pq_settings.write_json(o); });
        obj.entry("instance_file", instance_file);
    }
};

struct SharedData {
    KnapsackInstance<data_type> instance;
    std::atomic<data_type> solution{0};
    parallel_search::Termination termination;
};

struct Counter {
    long long pushed_nodes{0};
    long long processed_nodes{0};
    long long ignored_nodes{0};

    long long node_count() const noexcept {
        return pushed_nodes - processed_nodes - ignored_nodes;
    }
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

struct ThreadResult {
    Counter counter{};
    benchmark::Interval interval{};
};

[[gnu::noinline]] ThreadResult benchmark_thread(thread_coordination::Context& thread_context, pq_type& pq,
                                                SharedData& data) {
    ThreadResult result{};
    handle_type handle = pq.get_handle();
    if (thread_context.id() == 0) {
        auto [lb, ub] = data.instance.compute_bounds_linear(data.instance.capacity(), 0);
        data.solution.store(lb, std::memory_order_relaxed);
        handle.push(to_payload(ub, 0, data.instance.capacity(), 0));
        ++result.counter.pushed_nodes;
    }
    result.interval = data.termination.run(
        thread_context, handle, [&](node_type const& node) { process_node(node, handle, result.counter, data); },
        [&result]() { return result.counter.node_count(); });
    return result;
}

void run_benchmark(Settings const& settings) {
    std::clog << "Reading instance...\n";
    KnapsackInstance<data_type> instance;
    try {
        instance = KnapsackInstance<data_type>(settings.instance_file);
    } catch (std::runtime_error const& e) {
        std::cerr << "Error: " << settings.instance_file.string() << ": " << e.what() << '\n';
        std::exit(EXIT_FAILURE);
    }
    std::clog << "Instance has " << instance.size() << " items and " << std::fixed << instance.capacity()
              << " capacity\n";
    SharedData shared_data{std::move(instance), 0, parallel_search::Termination{settings.base_settings.num_threads}};
    std::vector<Counter> thread_counter(static_cast<std::size_t>(settings.base_settings.num_threads));
    std::vector<benchmark::Interval> thread_interval(static_cast<std::size_t>(settings.base_settings.num_threads));
    auto pq = pq_type(settings.base_settings.num_threads, std::size_t(10'000'000), settings.pq_settings);
    std::clog << "Working...\n";
    auto memory_start = memory_stats::Snapshot::take();
    thread_coordination::dispatch(settings.base_settings.affinity, settings.base_settings.num_threads, [&](auto ctx) {
        auto t_id = static_cast<std::size_t>(ctx.id());
        auto r = benchmark_thread(ctx, pq, shared_data);
        thread_counter[t_id] = r.counter;
        thread_interval[t_id] = r.interval;
    });
    auto memory_end = memory_stats::Snapshot::take();
    std::clog << "Done\n";
    Counter summed{};
    for (auto const& c : thread_counter) {
        summed.pushed_nodes += c.pushed_nodes;
        summed.processed_nodes += c.processed_nodes;
        summed.ignored_nodes += c.ignored_nodes;
    }
    assert(summed.node_count() == 0);
    std::clog << '\n';
    std::clog << "= Results =\n";
    std::clog << "Time (s): " << std::fixed << std::setprecision(3) << benchmark::seconds(thread_interval) << '\n';
    std::clog << "Solution: " << shared_data.solution.load() << '\n';
    std::clog << "Processed nodes: " << summed.processed_nodes << '\n';
    std::clog << "Ignored nodes: " << summed.ignored_nodes << '\n';
    {
        json::Object root{std::cout};
        root.object("settings", [&settings](json::Object& obj) { settings.write_json(obj); });
        root.object("instance", [&shared_data](json::Object& obj) {
            obj.entry("num_items", shared_data.instance.size());
            obj.entry("capacity", shared_data.instance.capacity());
        });
        root.object("results", [&](json::Object& results) {
            benchmark::write_timing(results, "", thread_interval);
            results.object("memory", [&](json::Object& memory) {
                memory_stats::write_json(memory, {{"start", memory_start}, {"end", memory_end}});
            });
            results.entry("processed_nodes", summed.processed_nodes);
            results.entry("ignored_nodes", summed.ignored_nodes);
            results.entry("solution", shared_data.solution.load());
        });
    }
    std::cout << '\n';
}

int main(int argc, char* argv[]) {
    return benchmark::run<pq_type, Settings>(argc, argv, run_benchmark);
}
