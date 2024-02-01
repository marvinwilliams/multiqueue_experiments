#include "util/build_info.hpp"
#include "util/knapsack_instance.hpp"

#include "cxxopts.hpp"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <queue>
#include <vector>

#ifdef FLOAT_INSTANCE
using data_type = double;
#else
using data_type = unsigned long;
#endif

struct Node {
    data_type upper_bound;
    std::size_t index;
    data_type free_capacity;
    data_type value;
};

bool operator<(Node const& lhs, Node const& rhs) noexcept {
    return lhs.upper_bound < rhs.upper_bound;
}

struct Settings {
    std::filesystem::path instance_file;
};

void register_cmd_options(Settings& settings, cxxopts::Options& cmd) {
    cmd.add_options()("instance", "The instance file", cxxopts::value<std::filesystem::path>(settings.instance_file),
                      "PATH");
    cmd.parse_positional({"instance"});
}

void write_settings_human_readable(Settings const& settings, std::ostream& out) {
    out << "Instance file: " << settings.instance_file << '\n';
}

void write_settings_json(Settings const& settings, std::ostream& out) {
    out << '{';
    out << std::quoted("instance_file") << ':' << settings.instance_file;
    out << '}';
}

void knapsack(Settings const& settings) noexcept {
    data_type best_value{0};
    long long processed_nodes{0};
    KnapsackInstance<data_type> instance;
    std::clog << "Reading instance...\n";
    try {
        instance = KnapsackInstance<data_type>(settings.instance_file);
    } catch (std::exception const& e) {
        std::cerr << "Error reading instance file: " << e.what() << '\n';
        std::exit(EXIT_FAILURE);
    }
    std::clog << "Instance has " << instance.size() << " items and " << std::fixed << instance.capacity() << " capacity\n";
    std::vector<Node> container;
    container.reserve(1 << 24);
    std::priority_queue<Node, std::vector<Node>> pq({}, std::move(container));
    std::clog << "Working...\n";
    auto t_start = std::chrono::steady_clock::now();
    {
        auto [lb, ub] = instance.compute_bounds_linear(instance.capacity(), 0);
        best_value = lb;
        pq.push(Node{ub, 0, instance.capacity(), 0});
    }
    while (!pq.empty()) {
        auto node = pq.top();
        pq.pop();
        ++processed_nodes;
        if (node.upper_bound <= best_value) {
            break;
        }
        auto [lb, ub] = instance.compute_bounds_linear(node.free_capacity, node.index + 1);
        if (node.value + lb > best_value) {
            best_value = node.value + lb;
        }
        if (node.index + 2 < instance.size()) {
            if (node.value + ub > best_value) {
                pq.push({node.value + ub, node.index + 1, node.free_capacity, node.value});
            }
            if (node.free_capacity >= instance.weight(node.index)) {
                node.value += instance.value(node.index);
                node.free_capacity -= instance.weight(node.index);
                ++node.index;
                pq.push(node);
            }
        }
    }
    auto t_end = std::chrono::steady_clock::now();
    std::clog << "Done\n\n";
    std::clog << "= Results =\n";
    std::clog << "Time (s): " << std::fixed << std::setprecision(3)
              << std::chrono::duration<double>(t_end - t_start).count() << '\n';
    std::clog << "Solution: " << best_value << '\n';
    std::clog << "Processed nodes: " << processed_nodes << '\n';

    std::cout << '{';
    std::cout << std::quoted("settings") << ':';
    write_settings_json(settings, std::cout);
    std::cout << ',';
    std::cout << std::quoted("instance") << ':';
    std::cout << '{';
    std::cout << std::quoted("num_items") << ':' << instance.size() << ',';
    std::cout << std::quoted("capacity") << ':' << std::fixed << instance.capacity();
    std::cout << '}' << ',';
    std::cout << std::quoted("results") << ':';
    std::cout << '{';
    std::cout << std::quoted("time_ns") << ':' << std::chrono::nanoseconds{t_end - t_start}.count() << ',';
    std::cout << std::quoted("processed_nodes") << ':' << processed_nodes << ',';
    std::cout << std::quoted("solution") << ':' << best_value;
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
    knapsack(settings);
    return EXIT_SUCCESS;
}
