#include "util/benchmark.hpp"
#include "util/graph.hpp"
#include "util/memory_stats.hpp"
#include "util/parallel_search.hpp"
#include "util/thread_coordination.hpp"
#include "wrapper/util/selector.hpp"

#include <cxxopts.hpp>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <atomic>
#include <cassert>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using pq_type = PQ<true, unsigned long, unsigned long>;
using handle_type = pq_type::handle_type;
using node_type = pq_type::value_type;

struct Settings {
    benchmark::BaseSettings base_settings{};
    pq_type::settings_type pq_settings{};
    std::filesystem::path graph_file;
    std::string source = "0";

    void register_cmd_options(cxxopts::Options& cmd) {
        base_settings.register_cmd_options(cmd);
        // clang-format off
        cmd.add_options()
            ("graph", "The input graph", cxxopts::value<std::filesystem::path>(graph_file), "PATH")
            ("source", "Source node: 0-based index or max-degree", cxxopts::value<std::string>(source), "NODE");
        // clang-format on
        pq_settings.register_cmd_options(cmd);
        cmd.parse_positional({"graph"});
    }

    bool validate() const {
        if (!base_settings.validate()) {
            return false;
        }
        if (!pq_settings.validate()) {
            return false;
        }
        if (graph_file.empty()) {
            std::cerr << "Error: No graph file specified\n";
            return false;
        }
        return true;
    }
    void write_human_readable(std::ostream& out) const {
        base_settings.write_human_readable(out);
        pq_settings.write_human_readable(out);
        out << "Graph: " << graph_file << '\n';
        out << "Source: " << source << '\n';
    }

    void write_json(json::Object& obj) const {
        base_settings.write_json(obj);
        obj.object("pq", [this](json::Object& o) { pq_settings.write_json(o); });
        obj.entry("graph", graph_file);
        obj.entry("source", source);
    }
};

struct alignas(L1_CACHE_LINE_SIZE) AtomicDistance {
    std::atomic<long long> value{std::numeric_limits<long long>::max()};
};

struct SharedData {
    Graph graph;
    std::vector<AtomicDistance> distances;
    std::size_t source{};
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
    auto current_distance = data.distances[node.second].value.load(std::memory_order_relaxed);
    if (static_cast<long long>(node.first) > current_distance) {
        ++counter.ignored_nodes;
        return;
    }
    for (auto i = data.graph.nodes[node.second]; i < data.graph.nodes[node.second + 1]; ++i) {
        auto target = data.graph.edges[i].target;
        auto d = static_cast<long long>(node.first) + data.graph.edges[i].weight;
        auto old_d = data.distances[target].value.load(std::memory_order_relaxed);
        while (d < old_d) {
            if (data.distances[target].value.compare_exchange_weak(old_d, d, std::memory_order_relaxed)) {
                if (handle.push({d, target})) {
                    ++counter.pushed_nodes;
                }
                break;
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
    auto handle = pq.get_handle();
    if (thread_context.id() == 0) {
        data.distances[data.source].value = 0;
        handle.push({0, data.source});
        ++result.counter.pushed_nodes;
    }
    result.interval = data.termination.run(
        thread_context, handle, [&](node_type const& node) { process_node(node, handle, result.counter, data); },
        [&result]() { return result.counter.node_count(); });
    return result;
}

void run_benchmark(Settings const& settings) {
    std::clog << "Reading graph...\n";
    SharedData shared_data{{}, {}, 0, parallel_search::Termination{settings.base_settings.num_threads}};
    try {
        shared_data.graph = Graph(settings.graph_file);
    } catch (std::runtime_error const& e) {
        std::cerr << "Error: " << settings.graph_file.string() << ": " << e.what() << '\n';
        std::exit(EXIT_FAILURE);
    }
    if (shared_data.graph.num_nodes() == 0) {
        std::cerr << "Error: " << settings.graph_file.string() << ": Graph has no nodes\n";
        std::exit(EXIT_FAILURE);
    }
    try {
        shared_data.source = shared_data.graph.source(settings.source);
    } catch (std::logic_error const& e) {
        std::cerr << "Error: --source " << settings.source << ": " << e.what() << '\n';
        std::exit(EXIT_FAILURE);
    }
    std::clog << "Graph has " << shared_data.graph.num_nodes() << " nodes and " << shared_data.graph.num_edges()
              << " edges, source " << shared_data.source << " has degree " << shared_data.graph.degree(shared_data.source)
              << '\n';
    shared_data.distances = std::vector<AtomicDistance>(shared_data.graph.num_nodes());

    std::vector<Counter> thread_counter(static_cast<std::size_t>(settings.base_settings.num_threads));
    std::vector<benchmark::Interval> thread_interval(static_cast<std::size_t>(settings.base_settings.num_threads));
    auto pq = pq_type(settings.base_settings.num_threads, shared_data.graph.num_nodes(), settings.pq_settings);
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
    auto furthest_node =
        std::max_element(shared_data.distances.begin(), shared_data.distances.end(), [](auto const& a, auto const& b) {
            auto a_val = a.value.load(std::memory_order_relaxed);
            auto b_val = b.value.load(std::memory_order_relaxed);
            if (b_val == std::numeric_limits<long long>::max()) {
                return false;
            }
            if (a_val == std::numeric_limits<long long>::max()) {
                return true;
            }
            return a_val < b_val;
        });
    std::clog << "= Results =\n";
    std::clog << "Time (s): " << std::fixed << std::setprecision(3) << benchmark::seconds(thread_interval) << '\n';
    std::clog << "Furthest node: " << furthest_node - shared_data.distances.begin() << '\n';
    std::clog << "Longest distance: " << furthest_node->value.load(std::memory_order_relaxed) << '\n';
    std::clog << "Processed nodes: " << summed.processed_nodes << '\n';
    std::clog << "Ignored nodes: " << summed.ignored_nodes << '\n';
    {
        json::Object root{std::cout};
        root.object("settings", [&settings](json::Object& obj) { settings.write_json(obj); });
        root.object("graph", [&shared_data](json::Object& graph) {
            graph.entry("num_nodes", shared_data.graph.num_nodes());
            graph.entry("num_edges", shared_data.graph.num_edges());
            graph.entry("source", shared_data.source);
            graph.entry("source_degree", shared_data.graph.degree(shared_data.source));
        });
        root.object("results", [&](json::Object& results) {
            results.entry("time_ns", benchmark::time_ns(thread_interval));
            results.entry("pushed_nodes", summed.pushed_nodes);
            results.object("memory", [&](json::Object& memory) {
                memory_stats::write_json(memory, {{"start", memory_start}, {"end", memory_end}});
            });
            results.entry("furthest_node", furthest_node - shared_data.distances.begin());
            results.entry("longest_distance", furthest_node->value.load(std::memory_order_relaxed));
            results.entry("processed_nodes", summed.processed_nodes);
            results.entry("ignored_nodes", summed.ignored_nodes);
            auto origin = benchmark::span(thread_interval).start;
            results.array("thread_data", thread_counter.begin(), thread_counter.end(),
                          [&, t = std::size_t{0}](std::ostream& out, Counter const& counter) mutable {
                              json::Object obj{out};
                              benchmark::write_interval(obj, "", thread_interval[t++], origin);
                              obj.entry("pushed_nodes", counter.pushed_nodes);
                              obj.entry("processed_nodes", counter.processed_nodes);
                              obj.entry("ignored_nodes", counter.ignored_nodes);
                          });
        });
    }
    std::cout << '\n';
}

int main(int argc, char* argv[]) {
    return benchmark::run<pq_type, Settings>(argc, argv, run_benchmark);
}
