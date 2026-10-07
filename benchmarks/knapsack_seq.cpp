#include "util/build_info.hpp"
#include "util/json.hpp"
#include "util/knapsack_instance.hpp"

#include "cxxopts.hpp"

#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
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
    void register_cmd_options(cxxopts::Options& cmd) {
        cmd.add_options()("instance", "The instance file", cxxopts::value<std::filesystem::path>(instance_file),
                          "PATH");
        cmd.parse_positional({"instance"});
    }

    bool validate() const {
        if (instance_file.empty()) {
            std::cerr << "Error: No instance file specified\n";
            return false;
        }
        return true;
    }
    void write_human_readable(std::ostream& out) const {
        out << "Instance file: " << instance_file << '\n';
    }

    void write_json(json::Object& obj) const {
        obj.entry("instance_file", instance_file);
    }
};

void knapsack(Settings const& settings) noexcept {
    data_type best_value{0};
    long long processed_nodes{0};
    std::size_t sum_sizes{0};
    std::size_t max_size{0};
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
        sum_sizes += pq.size();
        pq.pop();
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
        max_size = std::max(max_size, pq.size());
        ++processed_nodes;
    }
    auto t_end = std::chrono::steady_clock::now();
    std::clog << "Done\n\n";
    std::clog << "= Results =\n";
    std::clog << "Time (s): " << std::fixed << std::setprecision(3)
              << std::chrono::duration<double>(t_end - t_start).count() << '\n';
    std::clog << "Solution: " << best_value << '\n';
    std::clog << "Processed nodes: " << processed_nodes << '\n';
    auto average_pq_size =
        processed_nodes == 0 ? 0.0 : static_cast<double>(sum_sizes) / static_cast<double>(processed_nodes);
    std::clog << "Average PQ size: " << average_pq_size << '\n';
    std::clog << "Max PQ size: " << max_size << '\n';

    {
        json::Object root{std::cout};
        root.object("settings", [&settings](json::Object& obj) { settings.write_json(obj); });
        root.object("instance", [&instance](json::Object& obj) {
            obj.entry("num_items", instance.size());
            obj.entry("capacity", instance.capacity());
        });
        root.object("results", [&](json::Object& results) {
            results.array("thread_start_offset_ns", std::vector<long long>{0});
            results.array("thread_end_offset_ns", std::vector<long long>{std::chrono::nanoseconds{t_end - t_start}.count()});
            results.entry("processed_nodes", processed_nodes);
            results.entry("solution", best_value);
            results.entry("average_pq_size", average_pq_size);
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
    cmd.add_options()("h,help", "Print this help");
    Settings settings{};
    settings.register_cmd_options(cmd);

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
    settings.write_human_readable(std::clog);
    std::clog << '\n';
    if (!settings.validate()) {
        return EXIT_FAILURE;
    }
    std::clog << "= Running benchmark =\n";
    knapsack(settings);
    return EXIT_SUCCESS;
}
