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
static constexpr data_type default_min_weight = 0.01;
static constexpr data_type default_max_weight = 1.01;
static constexpr data_type default_min_add = 0.1;
static constexpr data_type default_max_add = 0.125;
#else
using data_type = unsigned long;
static constexpr data_type default_min_weight = 1000;
static constexpr data_type default_max_weight = 100000;
static constexpr data_type default_min_add = 10000;
static constexpr data_type default_max_add = 12500;
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
    long long n = 1000;
    data_type min_weight = default_min_weight;
    data_type max_weight = default_max_weight;
    data_type min_add = default_min_add;
    data_type max_add = default_max_add;
    double capacity_factor = 0.5;
    unsigned int seed = 1;
    std::filesystem::path instance_file;
};

void register_cmd_options(Settings& settings, cxxopts::Options& cmd) {
    // clang-format off
    cmd.add_options()
        ("n,num-elements", "Number of elements", cxxopts::value<long long>(settings.n), "NUMBER")
        ("a,min-weight", "Min weight", cxxopts::value<data_type>(settings.min_weight), "NUMBER")
        ("b,max-weight", "Max weight", cxxopts::value<data_type>(settings.max_weight), "NUMBER")
        ("l,min-add", "Min add to profits", cxxopts::value<data_type>(settings.min_add), "NUMBER")
        ("u,max-add", "Max add to profits", cxxopts::value<data_type>(settings.max_add), "NUMBER")
        ("f,factor", "Capacity as factor of expected total weight", cxxopts::value<double>(settings.capacity_factor), "NUMBER")
        ("s,seed", "Seed", cxxopts::value<unsigned int>(settings.seed), "NUMBER")
        ("instance", "The instance file", cxxopts::value<std::filesystem::path>(settings.instance_file), "PATH");
    // clang-format on
    cmd.parse_positional({"instance"});
}

void write_settings_human_readable(Settings const& settings, std::ostream& out) {
    if (settings.instance_file.empty()) {
        out << "Number of items: " << settings.n << '\n';
        out << "Weights: [" << settings.min_weight << ", " << settings.max_weight << "]\n";
        out << "Add to profit: [" << settings.min_add << ", " << settings.max_add << "]\n";
        out << "Capacity factor: " << settings.capacity_factor << '\n';
        out << "Seed: " << settings.seed << '\n';
    } else {
        out << "Instance file: " << settings.instance_file << '\n';
    }
}

void write_settings_json(Settings const& settings, std::ostream& out) {
    out << '{';
    if (settings.instance_file.empty()) {
        out << std::quoted("instance_type") << ':' << std::quoted("generated") << ',';
        out << std::quoted("num_elements") << ':' << settings.n << ',';
        out << std::quoted("min_weight") << ':' << settings.min_weight << ',';
        out << std::quoted("max_weight") << ':' << settings.max_weight << ',';
        out << std::quoted("min_add") << ':' << settings.min_add << ',';
        out << std::quoted("max_add") << ':' << settings.max_add << ',';
        out << std::quoted("capacity_factor") << ':' << settings.capacity_factor << ',';
        out << std::quoted("seed") << ':' << settings.seed;
    } else {
        out << std::quoted("instance_type") << ':' << std::quoted("file") << ',';
        out << std::quoted("instance_file") << ':' << settings.instance_file;
    }
    out << '}';
}

void knapsack(Settings const& settings) noexcept {
    data_type best_value{0};
    long long processed_nodes{0};
    KnapsackInstance<data_type> instance;
    if (settings.instance_file.empty()) {
        std::clog << "Generating instance...\n";
        instance = KnapsackInstance<data_type>(settings.n, settings.min_weight, settings.max_weight, settings.min_add,
                                               settings.max_add, settings.capacity_factor, settings.seed);
    } else {
        std::clog << "Reading instance...\n";
        try {
            instance = KnapsackInstance<data_type>(settings.instance_file);
        } catch (std::exception const& e) {
            std::cerr << "Error reading instance file: " << e.what() << '\n';
            std::exit(EXIT_FAILURE);
        }
    }
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

    std::clog << "= Running benchmark =\n";
    knapsack(settings);
    return EXIT_SUCCESS;
}
