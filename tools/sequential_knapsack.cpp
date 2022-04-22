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

using weight_type = unsigned long;
using value_type = unsigned long;

struct KnapsackInstance {
    struct Item {
        weight_type weight;
        value_type value;
    };
    std::vector<Item> items;
    weight_type capacity;
};

struct Node {
    value_type upper_bound;
    unsigned int index;
    weight_type free_capacity;
    value_type value;
};

bool operator<(Node const& lhs, Node const& rhs) noexcept {
    return lhs.upper_bound < rhs.upper_bound;
}

struct StatCounters {
    std::size_t pushed_nodes = 0;
    std::size_t ignored_nodes = 0;
    std::size_t extracted_nodes = 0;
    std::size_t processed_nodes = 0;
};

KnapsackInstance instance;
std::priority_queue<Node> pq;
value_type best_value;
StatCounters stats;

value_type get_lower_bound(weight_type capacity, unsigned int index) noexcept {
    value_type value = 0;
    while (index < instance.items.size() &&
           instance.items[index].weight <= capacity) {
        capacity -= instance.items[index].weight;
        value += instance.items[index].value;
        ++index;
    }
    return value;
}

value_type get_upper_bound(weight_type capacity, std::size_t index) noexcept {
    value_type result = 0;
    while (index < instance.items.size() &&
           capacity >= instance.items[index].weight) {
        result += instance.items[index].value;
        capacity -= instance.items[index].weight;
        ++index;
    }
    if (index < instance.items.size()) {
        result += static_cast<value_type>(
            static_cast<double>(capacity * instance.items[index].value) /
            static_cast<double>(instance.items[index].weight));
    }
    return result;
}

void read_problem(std::filesystem::path instance_file) {
    std::ifstream file_stream{instance_file};
    if (!file_stream) {
        throw std::runtime_error{"Could not open knapsack file"};
    }
    size_t n;
    file_stream >> n;
    if (!file_stream || file_stream.eof()) {
        throw std::runtime_error{"Error reading knapsack file"};
    }
    instance.items.clear();
    instance.items.reserve(n);
    file_stream >> instance.capacity;
    if (!file_stream || (n > 0 && file_stream.eof())) {
        throw std::runtime_error{"Error reading knapsack file"};
    }

    for (size_t i = 0; i < n; ++i) {
        if (!file_stream || file_stream.eof()) {
            throw std::runtime_error{"Error reading knapsack file"};
        }
        auto item = KnapsackInstance::Item{};
        file_stream >> item.value;
        if (!file_stream || file_stream.eof()) {
            throw std::runtime_error{"Error reading knapsack file"};
        }
        file_stream >> item.weight;
        if (!file_stream) {
            throw std::runtime_error{"Error reading knapsack file"};
        }
        instance.items.push_back(item);
    }
    std::sort(instance.items.begin(), instance.items.end(),
              [](auto const& lhs, auto const& rhs) {
                  return (static_cast<double>(lhs.value) /
                          static_cast<double>(lhs.weight)) >
                         (static_cast<double>(rhs.value) /
                          static_cast<double>(rhs.weight));
              });
}

void process_node(Node const& node) noexcept {
    if (node.index == instance.items.size()) {
        return;
    }
    if (node.upper_bound <= best_value) {
        // The upper bound of this node is worse than the currently best value
        ++stats.ignored_nodes;
        return;
    }
    ++stats.processed_nodes;

    // Check if there is enough capacity for the next item
    if (instance.items[node.index].weight <= node.free_capacity) {
        value_type new_value = node.value + instance.items[node.index].value;
        value_type new_capacity =
            node.free_capacity - instance.items[node.index].weight;
        if (new_value > best_value) {
            best_value = new_value;
        }
        pq.push({node.upper_bound, node.index + 1, new_capacity, new_value});
        ++stats.pushed_nodes;
    }
    value_type upper_bound_without_next =
        node.value + get_upper_bound(node.free_capacity, node.index + 1);
    if (upper_bound_without_next <= best_value) {
        return;
    }
    pq.push({upper_bound_without_next, node.index + 1, node.free_capacity,
             node.value});
    ++stats.pushed_nodes;
}

void main_loop() noexcept {
    while (!pq.empty()) {
        auto node = pq.top();
        pq.pop();
        ++stats.extracted_nodes;
        process_node(node);
    }
}

int main(int argc, char* argv[]) {
    std::filesystem::path instance_file;
    std::clog << "Command line: ";
    std::copy(argv, argv + argc,
              std::ostream_iterator<char const*>(std::clog, " "));
    std::clog << "\n\n";

    cxxopts::Options options("Sequential knapsack", "");
    // clang-format off
    options.add_options()
      ("f,file", "The input graph", cxxopts::value<std::filesystem::path>(instance_file)->default_value("knapsack.kp"), "PATH")
      ("h,help", "Print this help");
    // clang-format on

    try {
        auto result = options.parse(argc, argv);
        if (result.count("help")) {
            std::cerr << options.help() << std::endl;
            return 0;
        }
    } catch (cxxopts::OptionParseException const& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    std::clog << "Settings: \n\t"
              << "Instance file: " << instance_file.string();
    std::clog << "\n\n";

    std::clog << "Reading problem..." << std::flush;
    try {
        read_problem(instance_file);
    } catch (std::runtime_error const& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
    std::clog << "done\n\n";

    std::cout << "items: " << instance.items.size() << '\n';
    std::cout << "capacity: " << instance.capacity << '\n';
    value_type upper_bound = get_upper_bound(instance.capacity, 0);
    best_value = get_lower_bound(instance.capacity, 0);
    std::clog << "Solving knapsack instance..." << std::flush;
    if (best_value < upper_bound) {
        process_node({upper_bound, 0, instance.capacity, 0});
    }
    auto start = std::chrono::steady_clock::now();
    main_loop();
    auto end = std::chrono::steady_clock::now();
    std::clog << "done\n";
    std::cout << "value: " << best_value << '\n';
    std::cout << "time: " << std::setprecision(3)
              << std::chrono::duration<double>(end - start).count() << '\n';
    std::cout << "pushed nodes: " << stats.pushed_nodes << '\n';
    std::cout << "ignored nodes: " << stats.ignored_nodes << '\n';
    std::cout << "extracted nodes: " << stats.extracted_nodes << '\n';
    std::cout << "processed nodes: " << stats.processed_nodes << '\n';
}
