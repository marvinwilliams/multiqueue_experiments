#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "cxxopts.hpp"

struct Graph {
    struct Edge {
        std::size_t target;
        std::uint64_t weight;
    };
    std::vector<std::size_t> nodes;
    std::vector<Edge> edges;
};

struct Node {
    std::uint64_t distance;
    std::size_t id;
};

bool operator>(Node const& lhs, Node const& rhs) noexcept {
    return lhs.distance > rhs.distance;
}

struct StatCounters {
    std::size_t pushed_nodes = 0;
    std::size_t ignored_nodes = 0;
    std::size_t extracted_nodes = 0;
    std::size_t processed_nodes = 0;
};

Graph graph;
std::priority_queue<Node, std::vector<Node>, std::greater<>> pq;
std::unique_ptr<std::uint64_t[]> shortest_distances;
StatCounters stats;

void read_graph(std::filesystem::path const& graph_file) {
    std::ifstream graph_stream{graph_file};
    if (!graph_stream) {
        throw std::runtime_error{"Could not open graph file"};
    }
    std::vector<std::vector<Graph::Edge>> edges_per_node;
    std::size_t source;
    std::size_t target;
    std::uint32_t weight;
    std::string problem;
    std::string first;
    while (graph_stream >> first) {
        if (first == "c") {
            graph_stream.ignore(std::numeric_limits<std::streamsize>::max(),
                                '\n');
        } else if (first == "p") {
            graph_stream >> problem;
            std::size_t num_nodes;
            std::size_t num_edges;
            graph_stream >> num_nodes >> num_edges;
            graph.nodes.resize(num_nodes + 1, 0);
            edges_per_node.resize(num_nodes);
            graph.edges.reserve(num_edges);
        } else if (first == "a") {
            graph_stream >> source >> target >> weight;
            edges_per_node[source - 1].push_back(
                Graph::Edge{target - 1, weight});
        } else {
            throw std::runtime_error{"Error reading file"};
        }
    }
    for (std::size_t i = 0; i < edges_per_node.size(); ++i) {
        graph.nodes[i + 1] = graph.nodes[i] + edges_per_node[i].size();
        std::copy(edges_per_node[i].begin(), edges_per_node[i].end(),
                  std::back_inserter(graph.edges));
    }
}

void main_loop() noexcept {
    while (!pq.empty()) {
        auto node = pq.top();
        pq.pop();
        ++stats.extracted_nodes;
        // Ignore stale nodes
        if (node.distance > shortest_distances[node.id]) {
            ++stats.ignored_nodes;
            continue;
        }
        ++stats.processed_nodes;
        for (std::size_t i = graph.nodes[node.id]; i < graph.nodes[node.id + 1];
             ++i) {
            std::uint64_t d = node.distance + graph.edges[i].weight;
            if (d < shortest_distances[graph.edges[i].target]) {
                shortest_distances[graph.edges[i].target] = d;
                pq.push({d, graph.edges[i].target});
                ++stats.pushed_nodes;
            }
        }
    }
}

int main(int argc, char* argv[]) {
    std::clog << "Command line: ";
    std::copy(argv, argv + argc,
              std::ostream_iterator<char const*>(std::clog, " "));

    std::clog << "\n\n";

    std::filesystem::path graph_file;
    std::filesystem::path output;

    cxxopts::Options options("Sequential shortest path",
                             "Compute shortest paths");
    // clang-format off
    options.add_options()
      ("f,file", "The input graph", cxxopts::value<std::filesystem::path>(graph_file)->default_value("graph.gr"), "PATH")
      ("n,start", "The starting node"
       "(default: 0)", cxxopts::value<std::size_t>(), "NUMBER")
      ("o,output", "File to output the solution to (default: no output)", cxxopts::value<std::filesystem::path>(output), "PATH")
      ("h,help", "Print this help");
    // clang-format on

    std::size_t starting_node = 0;

    try {
        auto result = options.parse(argc, argv);
        if (result.count("help")) {
            std::cerr << options.help() << std::endl;
            return 0;
        }
        if (result.count("start") > 0) {
            starting_node = result["start"].as<std::size_t>();
        }
    } catch (cxxopts::OptionParseException const& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    std::clog << "Reading graph..." << std::flush;
    try {
        read_graph(graph_file);
    } catch (std::runtime_error const& e) {
        std::cerr << '\n' << e.what() << '\n';
        return 1;
    }

    std::clog << "done\n";
    std::cout << "nodes: " << graph.nodes.size() - 1 << '\n';
    std::cout << "edges: " << graph.edges.size() << '\n';

    shortest_distances =
        std::make_unique<std::uint64_t[]>(graph.nodes.size() - 1);
    for (std::size_t i = 0; i + 1 < graph.nodes.size(); ++i) {
        shortest_distances[i] = std::numeric_limits<std::uint64_t>::max();
    }
    std::clog << "\nComputing shortest paths..." << std::flush;
    shortest_distances[starting_node] = 0;
    pq.push({0, starting_node});
    auto start = std::chrono::steady_clock::now();
    main_loop();
    auto end = std::chrono::steady_clock::now();
    std::clog << "done\n";
    std::cout << "time: " << std::setprecision(3)
              << std::chrono::duration<double>(end - start).count() << '\n';
    std::cout << "pushed nodes: " << stats.pushed_nodes << '\n';
    std::cout << "ignored nodes: " << stats.ignored_nodes << '\n';
    std::cout << "extracted nodes: " << stats.extracted_nodes << '\n';
    std::cout << "processed nodes: " << stats.processed_nodes << '\n';
    if (!output.empty()) {
        std::clog << "\nWriting output..." << std::flush;
        std::ofstream out_stream{output};
        if (!out_stream) {
            std::cerr << "\nCould not open output file\n";
            return 1;
        }
        out_stream << "node dist\n";
        for (std::size_t i = 0; i + 1 < graph.nodes.size(); ++i) {
            out_stream << i << ' ' << shortest_distances[i] << '\n';
        }
        std::clog << "done\n";
    }
}
