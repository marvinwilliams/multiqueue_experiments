#include "build_info.hpp"
#include "graph.hpp"

#include "cxxopts.hpp"

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

using clock_type = std::chrono::steady_clock;

struct Node {
    long long distance;
    std::size_t id;

    friend bool operator>(Node const& lhs, Node const& rhs) noexcept {
        return lhs.distance > rhs.distance;
    }
    friend bool operator<(Node const& lhs, Node const& rhs) noexcept {
        return lhs.distance < rhs.distance;
    }
};

#ifdef USE_FIFO
using PriorityQueue = std::queue<Node>;
#else
#ifdef REVERSE_PRIORITY
using PriorityQueue = std::priority_queue<Node, std::vector<Node>, std::less<>>;
#else
using PriorityQueue = std::priority_queue<Node, std::vector<Node>, std::greater<>>;
#endif
#endif

struct Data {
    std::vector<long long> shortest_distances;
    long long processed_nodes{0};
    long long ignored_nodes{0};

    Data(std::size_t num_nodes) : shortest_distances(num_nodes, std::numeric_limits<long long>::max()) {
    }
};

void dijkstra(PriorityQueue& pq, Data& data, Graph const& graph) noexcept {
    while (!pq.empty()) {
#ifdef USE_FIFO
        auto node = pq.front();
#else
        auto node = pq.top();
#endif
        pq.pop();
        // Ignore stale nodes
        if (node.distance > data.shortest_distances[node.id]) {
            ++data.ignored_nodes;
            continue;
        }
        ++data.processed_nodes;
        for (std::size_t i = graph.nodes[node.id]; i < graph.nodes[node.id + 1]; ++i) {
            auto d = node.distance + graph.edges[i].weight;
            if (d < data.shortest_distances[graph.edges[i].target]) {
                data.shortest_distances[graph.edges[i].target] = d;
                pq.push({d, graph.edges[i].target});
            }
        }
    }
}

int main(int argc, char* argv[]) {
    write_build_info(std::clog);
    std::clog << "\nCommand line:";

    for (int i = 0; i < argc; ++i) {
        std::clog << ' ' << argv[i];
    }
    std::clog << '\n' << '\n';

    std::filesystem::path graph_file;
    std::filesystem::path distance_file;

    cxxopts::Options options(argv[0]);
    // clang-format off
    options.add_options()
      ("file", "The input graph", cxxopts::value<std::filesystem::path>(graph_file)->default_value("graph.gr"), "PATH")
      ("o,distance-file", "Path to write the distances to", cxxopts::value<std::filesystem::path>(distance_file), "PATH")
      ("h,help", "Print this help");
    // clang-format on
    options.parse_positional({"file"});

    {
        cxxopts::ParseResult result;
        try {
            result = options.parse(argc, argv);
        } catch (cxxopts::OptionParseException const& e) {
            std::cerr << "Error parsing arguments: " << e.what() << '\n';
            std::cerr << options.help() << std::endl;
            return 1;
        }
        if (result.count("help") > 0U) {
            std::cerr << options.help() << std::endl;
            return 0;
        }
    }
    std::ofstream distance_out;
    if (!distance_file.empty()) {
        distance_out = std::ofstream(distance_file);
        if (!distance_out) {
            std::cerr << "Error: Could not open file " << distance_file << " for writing" << std::endl;
            return EXIT_FAILURE;
        }
    }
    std::clog << "Reading graph..." << std::endl;
    Graph graph;
    try {
        graph = Graph(graph_file);
    } catch (std::runtime_error const& e) {
        std::cerr << "\nError reading graph: " << e.what() << std::endl;
        return 1;
    }

    PriorityQueue pq;
    Data data(graph.num_nodes());
    data.shortest_distances[0] = 0;
    pq.push({0, 0});

    std::clog << "Computing shortest paths..." << std::endl;
    auto t_start = std::chrono::steady_clock::now();
    dijkstra(pq, data, graph);
    auto t_end = std::chrono::steady_clock::now();
    auto time = std::chrono::duration<double>(t_end - t_start).count();

    if (distance_out.is_open()) {
        std::clog << "Writing distances..." << std::endl;
        for (std::size_t i = 0; i < graph.num_nodes(); ++i) {
            distance_out << i << ' ' << data.shortest_distances[i] << '\n';
        }
        distance_out.close();
    }
    std::clog << "Finished\n" << std::endl;

    std::clog << "Time (s): " << std::setprecision(3) << time << '\n';
    std::clog << "Processed nodes: " << data.processed_nodes << '\n';
    std::clog << "Ignored nodes: " << data.ignored_nodes << '\n';
    std::clog << "Total nodes: " << data.processed_nodes + data.ignored_nodes << '\n';

    std::cout << "graph,nodes,edges,time,processed_nodes,ignored_nodes\n";
    std::cout << graph_file.string() << ',' << graph.num_nodes() << ',' << graph.num_edges() << ',' << time << ','
              << data.processed_nodes << ',' << data.ignored_nodes << std::endl;
}
