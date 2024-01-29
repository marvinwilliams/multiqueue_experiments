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

using pq_type = std::priority_queue<Node>;

struct Settings {
    long long n = 1000;
    data_type min_weight = default_min_weight;
    data_type max_weight = default_max_weight;
    data_type min_add = default_min_add;
    data_type max_add = default_max_add;
    double capacity_factor = 0.5;
    unsigned int seed = 1;
};

void register_cmd_options(Settings& settings, cxxopts::Options& cmd) {
    cmd.add_options()
        // clang-format off
        ("n,num-elements", "Number of elements", cxxopts::value<long long>(settings.n), "NUMBER")
        ("a,min-weight", "Min weight", cxxopts::value<data_type>(settings.min_weight), "NUMBER")
        ("b,max-weight", "Max weight", cxxopts::value<data_type>(settings.max_weight), "NUMBER")
        ("l,min-add", "Min add to profits", cxxopts::value<data_type>(settings.min_add), "NUMBER")
        ("u,max-add", "Max add to profits", cxxopts::value<data_type>(settings.max_add), "NUMBER")
        ("f,factor", "Capacity as factor of expected total weight", cxxopts::value<double>(settings.capacity_factor), "NUMBER")
        ("s,seed", "Seed", cxxopts::value<unsigned int>(settings.seed), "NUMBER");
    // clang-format on
}

void write_settings_human_readable(Settings const& settings, std::ostream& out) {
    out << "Number of items: " << settings.n << '\n';
    out << "Weights: [" << settings.min_weight << ", " << settings.max_weight << "]\n";
    out << "Add to profit: [" << settings.min_add << ", " << settings.max_add << "]\n";
    out << "Capacity factor: " << settings.capacity_factor << '\n';
    out << "Seed: " << settings.seed << '\n';
}

void write_settings_json(Settings const& settings, std::ostream& out) {
    out << '{';
    out << std::quoted("num-elements") << ':' << settings.n << ',';
    out << std::quoted("min-weight") << ':' << settings.min_weight << ',';
    out << std::quoted("max-weight") << ':' << settings.max_weight << ',';
    out << std::quoted("min-add") << ':' << settings.min_add << ',';
    out << std::quoted("max-add") << ':' << settings.max_add << ',';
    out << std::quoted("capacity-factor") << ':' << settings.capacity_factor << ',';
    out << std::quoted("seed") << ':' << settings.seed;
    out << '}';
}

void knapsack(Settings const& settings) noexcept {
    data_type best_value{0};
    long long processed_nodes{0};
    std::clog << "Generating instance...\n";
    auto instance = KnapsackInstance<data_type>(settings.n, settings.min_weight, settings.max_weight, settings.min_add,
                                                settings.max_add, settings.capacity_factor, settings.seed);
    pq_type::container_type container;
    container.reserve(1 << 24);
    pq_type pq{{}, std::move(container)};
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
    std::clog << "Done\n" << std::endl;

    std::clog << "Time (s): " << std::fixed << std::setprecision(3)
              << std::chrono::duration<double>(t_end - t_start).count() << '\n';
    std::clog << "Solution: " << best_value << '\n';
    std::clog << "Processed nodes: " << processed_nodes << '\n';

    std::cout << '{';
    std::cout << std::quoted("settings") << ':';
    write_settings_json(settings, std::cout);
    std::cout << ',';
    std::cout << std::quoted("time-ns") << ':' << std::chrono::nanoseconds{t_end - t_start}.count() << ',';
    std::cout << std::quoted("processed-nodes") << ':' << processed_nodes << ',';
    std::cout << std::quoted("solution") << ':' << best_value;
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
            std::cerr << cmd.help() << std::endl;
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
