#include "knapsack_instance.hpp"

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

struct Node {
    value_type upper_bound;
    unsigned int index;
    weight_type free_capacity;
    value_type value;

    friend bool operator<(Node const& lhs, Node const& rhs) noexcept {
        return lhs.upper_bound < rhs.upper_bound;
    }
};

#ifdef USE_FIFO
using PriorityQueue = std::queue<Node>;
#else
using PriorityQueue = std::priority_queue<Node>;
#endif

struct Data {
    unsigned long long best_value{0};
    long long processed_nodes = 0;
    long long ignored_nodes = 0;
};

value_type lower_bound(KnapsackInstance const& instance, weight_type capacity, unsigned int index) noexcept {
    value_type value = 0;
    while (index < instance.items.size() && instance.items[index].weight <= capacity) {
        capacity -= instance.items[index].weight;
        value += instance.items[index].value;
        ++index;
    }
    return value;
}

unsigned long upper_bound(KnapsackInstance const& instance, unsigned long capacity, std::size_t index) noexcept {
    unsigned long value_offset = instance.prefix_sum[index].value;
    unsigned long target_capacity = instance.prefix_sum[index].weight + capacity;
    while (index != instance.prefix_sum.size()) {
        if (instance.prefix_sum[index].weight > target_capacity) {
            double fraction = static_cast<double>(target_capacity - instance.prefix_sum[index - 1].weight) /
                static_cast<double>(instance.items[index - 1].weight);
            return (instance.prefix_sum[index - 1].value - value_offset) +
                static_cast<unsigned long>(static_cast<double>(instance.items[index - 1].value) * fraction);
        }
        ++index;
    }
    // All items fit
    return instance.prefix_sum[index - 1].value - value_offset;
}

void knapsack(PriorityQueue& pq, Data& data, KnapsackInstance const& instance) noexcept {
    while (!pq.empty()) {
#ifdef USE_FIFO
        auto node = pq.front();
#else
        auto node = pq.top();
#endif
        pq.pop();
        if (node.upper_bound <= data.best_value) {
            // The upper bound of this node is worse than the currently best value
            /* std::cout << "Ignored " << node.upper_bound << '\n'; */
            ++data.ignored_nodes;
            continue;
        }
        /* std::cout << "Popped " << node.upper_bound << ' ' << node.index << ' ' << node.free_capacity << ' ' */
        /*           << node.value << '\n'; */
        ++data.processed_nodes;
        if (node.index == instance.items.size()) {
            if (node.value > data.best_value) {
                data.best_value = node.value;
            }
            continue;
        }

        // Check if there is enough capacity for the next item
        if (instance.items[node.index].weight <= node.free_capacity) {
            value_type new_value = node.value + instance.items[node.index].value;
            value_type new_capacity = node.free_capacity - instance.items[node.index].weight;
            pq.push({node.upper_bound, node.index + 1, new_capacity, new_value});
            /* std::cout << "Pushed " << node.upper_bound << ' ' << node.index + 1 << ' ' << new_capacity << ' ' */
            /*           << new_value << '\n'; */
        }
        value_type upper_bound_without_next = node.value + upper_bound(instance, node.free_capacity, node.index + 1);
        if (upper_bound_without_next <= data.best_value) {
            continue;
        }
        pq.push({upper_bound_without_next, node.index + 1, node.free_capacity, node.value});
        /* std::cout << "Pushed " << upper_bound_without_next << ' ' << node.index + 1 << ' ' << node.free_capacity << ' ' */
        /*           << node.value << '\n'; */
    }
}

void print_header() {
    std::clog << "Built on " << __DATE__ << ' ' << __TIME__ << " with:\n";
#ifdef NDEBUG
    std::clog << "  Release build\n";
#else
    std::clog << "  Debug build\n";
#endif
#ifdef __clang__
    std::clog << "  Clang " << __clang_version__ << '\n';
#elif defined(__GNUC__)
    std::clog << "  GCC " << __VERSION__ << '\n';
#else
    std::clog << "  Unknown compiler\n";
#endif
}

int main(int argc, char* argv[]) {
    print_header();
    std::clog << "\nCommand line:";

    for (int i = 0; i < argc; ++i) {
        std::clog << ' ' << argv[i];
    }
    std::clog << '\n' << '\n';

    std::filesystem::path instance_file;

    cxxopts::Options options(argv[0]);
    // clang-format off
    options.add_options()
      ("file", "The knapsack instance", cxxopts::value<std::filesystem::path>(instance_file), "PATH")
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

    std::clog << "Reading instance..." << std::flush;
    KnapsackInstance instance;
    try {
        instance.from_file(instance_file);
    } catch (std::runtime_error const& e) {
        std::clog << "failed: " << e.what() << std::endl;
        return 1;
    }

    std::clog << "done\n";

    PriorityQueue pq;
    Data data;

    auto ub = upper_bound(instance, instance.capacity, 0);
    auto lb = lower_bound(instance, instance.capacity, 0);
    /* std::cout << "Upper bound: " << ub << '\n'; */
    /* std::cout << "Lower bound: " << lb << '\n'; */
    data.best_value = lb;
    pq.push({ub, 0, instance.capacity, 0});
    /* std::cout << "Pushed " << ub << ' ' << 0 << ' ' << instance.capacity << ' ' << 0 << '\n'; */
    std::clog << "Solving knapsack instance..." << std::flush;
    auto t_start = std::chrono::steady_clock::now();
    knapsack(pq, data, instance);
    auto t_end = std::chrono::steady_clock::now();
    std::clog << "done\n";
    auto time = std::chrono::duration<double>(t_end - t_start).count();

    std::clog << "Time (s): " << std::setprecision(3) << time << '\n';
    std::clog << "Value: " << data.best_value << '\n';
    std::clog << "Processed nodes: " << data.processed_nodes << '\n';
    std::clog << "Ignored nodes: " << data.ignored_nodes << '\n';
    std::clog << "Total nodes: " << data.processed_nodes + data.ignored_nodes << '\n';

    std::cout << "instance,nodes,time,processed_nodes,ignored_nodes,value\n";
    std::cout << instance_file.string() << ',' << instance.items.size() << ',' << time << ',' << data.processed_nodes
              << ',' << data.ignored_nodes << ',' << data.best_value << std::endl;
}
