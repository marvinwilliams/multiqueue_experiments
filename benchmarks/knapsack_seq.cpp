#include "build_info.hpp"
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

using payload_type = unsigned long;

struct Node {
    long long upper_bound;
    std::size_t index;
    long long free_capacity;
    long long value;

    friend bool operator<(Node const& lhs, Node const& rhs) noexcept {
        return lhs.upper_bound < rhs.upper_bound;
    }
};

#ifdef USE_FIFO
using pq_type = std::queue<Node>;
#else
using pq_type = std::priority_queue<Node>;
#endif

struct Data {
    long long best_value{0};
    long long pushed_nodes{0};
    long long ignored_nodes{0};
    long long processed_nodes{0};
};

void knapsack(pq_type& pq, Data& data, KnapsackInstance const& instance) noexcept {
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
        if (node.index == instance.size()) {
            if (node.value > data.best_value) {
                data.best_value = node.value;
            }
            continue;
        }

        // Check if there is enough capacity for the next item
        if (instance.items()[node.index].weight <= node.free_capacity) {
            auto new_value = node.value + instance.items()[node.index].value;
            auto new_capacity = node.free_capacity - instance.items()[node.index].weight;
            pq.push({node.upper_bound, node.index + 1, new_capacity, new_value});
        }
        auto upper_bound_without_next = node.value + instance.upper_bound(node.free_capacity, node.index + 1);
        if (upper_bound_without_next <= data.best_value) {
            continue;
        }
        pq.push({upper_bound_without_next, node.index + 1, node.free_capacity, node.value});
    }
}

int main(int argc, char* argv[]) {
    write_build_info(std::clog);
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

    try {
        auto result = options.parse(argc, argv);
        if (result.count("help") > 0U) {
            std::cerr << options.help() << std::endl;
            return EXIT_SUCCESS;
        }
    } catch (cxxopts::OptionParseException const& e) {
        std::cerr << "Error parsing arguments: " << e.what() << '\n';
        std::cerr << options.help() << std::endl;
        return 1;
    }

    KnapsackInstance instance;
    try {
        instance = KnapsackInstance(instance_file);
    } catch (std::runtime_error const& e) {
        std::clog << "failed: " << e.what() << std::endl;
        return 1;
    }

    pq_type pq;
    Data data;

    auto ub = instance.upper_bound(instance.capacity(), 0);
    auto lb = instance.lower_bound(instance.capacity(), 0);
    data.best_value = lb;
    pq.push({ub, 0, instance.capacity(), 0});
    /* std::cout << "Pushed " << ub << ' ' << 0 << ' ' << instance.capacity << ' ' << 0 << '\n'; */
    std::clog << "Solving knapsack instance...\n";
    auto t_start = std::chrono::steady_clock::now();
    knapsack(pq, data, instance);
    auto t_end = std::chrono::steady_clock::now();
    std::clog << "Finished\n" << std::endl;

    std::clog << "Time (s): " << std::setprecision(3) << std::chrono::duration<double>(t_end - t_start).count() << '\n';
    std::clog << "Value: " << data.best_value << '\n';
    std::clog << "Processed nodes: " << data.processed_nodes << '\n';
    std::clog << "Ignored nodes: " << data.ignored_nodes << '\n';

    std::cout << "instance,nodes,time,processed_nodes,ignored_nodes,value\n";
    std::cout << instance_file.string() << ',' << instance.size() << ',' << (t_end - t_start).count() << ','
              << data.processed_nodes << ',' << data.ignored_nodes << ',' << data.best_value << std::endl;
}
