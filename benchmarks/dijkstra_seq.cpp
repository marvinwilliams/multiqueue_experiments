#include "util/build_info.hpp"
#include "util/graph.hpp"

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

void dijkstra(std::filesystem::path const& graph_file) noexcept {
    std::clog << "Reading graph...\n";
    Graph graph;
    try {
        graph = Graph(graph_file);
    } catch (std::runtime_error const& e) {
        std::clog << "Error: " << e.what() << '\n';
        std::exit(EXIT_FAILURE);
    }
    std::clog << "Graph has " << graph.num_nodes() << " nodes and " << graph.num_edges() << " edges\n";
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
    distances[0] = 0;
    pq.push({0, 0});
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
    auto longest_distance = *std::max_element(distances.begin(), distances.end(), [](auto const& a, auto const& b) {
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
    std::clog << "Longest distance: " << longest_distance << '\n';
    std::clog << "Processed nodes: " << processed_nodes << '\n';
    std::clog << "Ignored nodes: " << ignored_nodes << '\n';
    std::clog << "Average PQ size: " << static_cast<double>(sum_sizes) / static_cast<double>(processed_nodes + ignored_nodes) << '\n';
    std::clog << "Max PQ size: " << max_size << '\n';

    std::cout << '{';
    std::cout << std::quoted("settings") << ':';
    std::cout << '{';
    std::cout << std::quoted("graph_file") << ':' << graph_file;
    std::cout << '}' << ',';
    std::cout << std::quoted("graph") << ':';
    std::cout << '{';
    std::cout << std::quoted("num_nodes") << ':' << graph.num_nodes() << ',';
    std::cout << std::quoted("num_edges") << ':' << graph.num_edges();
    std::cout << '}' << ',';
    std::cout << std::quoted("results") << ':';
    std::cout << '{';
    std::cout << std::quoted("time_ns") << ':' << std::chrono::nanoseconds{t_end - t_start}.count() << ',';
    std::cout << std::quoted("longest_distance") << ':' << longest_distance << ',';
    std::cout << std::quoted("processed_nodes") << ':' << processed_nodes << ',';
    std::cout << std::quoted("ignored_nodes") << ':' << ignored_nodes << ',';
    std::cout << std::quoted("average_pq_size") << ':'
              << static_cast<double>(sum_sizes) / static_cast<double>(processed_nodes + ignored_nodes) << ',';
    std::cout << std::quoted("max_pq_size") << ':' << max_size;
    std::cout << '}';
    std::cout << '}' << '\n';
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
    // clang-format off
    cmd.add_options()
        ("h,help", "Print this help")
        ("graph", "The input graph", cxxopts::value<std::filesystem::path>(graph_file), "PATH");
    // clang-format on
    cmd.parse_positional({"graph"});

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
    std::clog << "Graph: " << graph_file << '\n';
    std::clog << '\n';

    std::clog << "= Running benchmark =\n";
    dijkstra(graph_file);
    return EXIT_SUCCESS;
}
