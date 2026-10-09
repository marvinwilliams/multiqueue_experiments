#include "util/build_info.hpp"
#include "util/graph.hpp"
#include "util/json.hpp"

#include <cxxopts.hpp>

#include <chrono>
#include <filesystem>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using clock_type = std::chrono::steady_clock;

struct Node {
    long long distance;
    std::size_t id;

    friend bool operator>(Node const& lhs, Node const& rhs) noexcept {
        return lhs.distance > rhs.distance;
    }
};

void dijkstra(std::filesystem::path const& graph_file, std::string const& source_spec) noexcept {
    std::clog << "Reading graph...\n";
    Graph graph;
    try {
        graph = Graph(graph_file);
    } catch (std::runtime_error const& e) {
        std::cerr << "Error: " << graph_file.string() << ": " << e.what() << '\n';
        std::exit(EXIT_FAILURE);
    }
    std::size_t source = 0;
    try {
        source = graph.source(source_spec);
    } catch (std::logic_error const& e) {
        std::cerr << "Error: --source " << source_spec << ": " << e.what() << '\n';
        std::exit(EXIT_FAILURE);
    }
    std::clog << "Graph has " << graph.num_nodes() << " nodes and " << graph.num_edges() << " edges, source " << source
              << " has degree " << graph.degree(source) << '\n';
    std::vector<long long> distances(graph.num_nodes(), std::numeric_limits<long long>::max());
    long long processed_nodes{0};
    long long ignored_nodes{0};
    std::size_t sum_sizes{0};
    std::size_t max_size{0};
    std::vector<Node> container;
    container.reserve(graph.num_nodes());
    std::priority_queue<Node, std::vector<Node>, std::greater<>> pq({}, std::move(container));
    std::clog << "Working...\n";
    auto t_start = std::chrono::steady_clock::now();
    distances[source] = 0;
    pq.push({0, source});
    while (!pq.empty()) {
        sum_sizes += pq.size();
        max_size = std::max(max_size, pq.size());
        auto node = pq.top();
        pq.pop();
        // Ignore stale nodes
        if (node.distance > distances[node.id]) {
            ++ignored_nodes;
            continue;
        }
        for (std::size_t i = graph.nodes[node.id]; i < graph.nodes[node.id + 1]; ++i) {
            auto d = node.distance + graph.edges[i].weight;
            if (d < distances[graph.edges[i].target]) {
                distances[graph.edges[i].target] = d;
                pq.push({d, graph.edges[i].target});
            }
        }
        ++processed_nodes;
    }
    auto t_end = std::chrono::steady_clock::now();
    std::clog << "Done\n\n";
    auto furthest_node = std::max_element(distances.begin(), distances.end(), [](auto const& a, auto const& b) {
        if (b == std::numeric_limits<long long>::max()) {
            return false;
        }
        if (a == std::numeric_limits<long long>::max()) {
            return true;
        }
        return a < b;
    });
    std::clog << "= Results =\n";
    std::clog << "Time (s): " << std::fixed << std::setprecision(3)
              << std::chrono::duration<double>(t_end - t_start).count() << '\n';
    std::clog << "Furthest node: " << furthest_node - distances.begin() << '\n';
    std::clog << "Longest distance: " << *furthest_node << '\n';
    std::clog << "Processed nodes: " << processed_nodes << '\n';
    std::clog << "Ignored nodes: " << ignored_nodes << '\n';
    std::clog << "Average PQ size: " << static_cast<double>(sum_sizes) / static_cast<double>(processed_nodes + ignored_nodes) << '\n';
    std::clog << "Max PQ size: " << max_size << '\n';

    {
        json::Object root{std::cout};
        root.object("settings", [&](json::Object& obj) {
            obj.entry("graph_file", graph_file);
            obj.entry("source", source_spec);
        });
        root.object("graph", [&](json::Object& obj) {
            obj.entry("num_nodes", graph.num_nodes());
            obj.entry("num_edges", graph.num_edges());
            obj.entry("source", source);
            obj.entry("source_degree", graph.degree(source));
        });
        root.object("results", [&](json::Object& results) {
            results.entry("time_ns", std::chrono::nanoseconds{t_end - t_start}.count());
            results.entry("furthest_node", furthest_node - distances.begin());
            results.entry("longest_distance", *furthest_node);
            results.entry("processed_nodes", processed_nodes);
            results.entry("ignored_nodes", ignored_nodes);
            results.entry("average_pq_size",
                          static_cast<double>(sum_sizes) / static_cast<double>(processed_nodes + ignored_nodes));
            results.entry("max_pq_size", max_size);
        });
    }
    std::cout << '\n';
}

int main(int argc, char* argv[]) {
    write_build_info(std::clog);
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
    std::filesystem::path graph_file;
    std::string source = "0";
    // clang-format off
    cmd.add_options()
        ("h,help", "Print this help")
        ("graph", "The input graph", cxxopts::value<std::filesystem::path>(graph_file), "PATH")
        ("source", "Source node: 0-based index or max-degree", cxxopts::value<std::string>(source), "NODE");
    // clang-format on
    cmd.parse_positional({"graph"});

    try {
        auto args = cmd.parse(argc, argv);
        if (args.count("help") > 0) {
            std::cerr << cmd.help() << '\n';
            return EXIT_SUCCESS;
        }
    } catch (cxxopts::OptionException const& e) {
        std::cerr << "Error parsing command line: " << e.what() << '\n';
        std::cerr << "Use --help for usage information" << '\n';
        return EXIT_FAILURE;
    }

    std::clog << "= Settings =\n";
    std::clog << "Graph: " << graph_file << '\n';
    std::clog << "Source: " << source << '\n';
    std::clog << '\n';

    std::clog << "= Running benchmark =\n";
    dijkstra(graph_file, source);
    return EXIT_SUCCESS;
}
