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
    /* std::size_t hint; */
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
        ++data.processed_nodes;
        if (node.upper_bound <= data.best_value) {
            return;
        }
        if (node.index + 1 == instance.size()) {
            continue;
        }
        /* auto hint = node.hint; */
        /* auto const& [lb, ub] = instance.compute_bounds_hint(node.free_capacity, node.index + 1, hint); */
        auto const& [lb, ub] = instance.compute_bounds_linear(node.free_capacity, node.index + 1);
        if (node.value + lb > data.best_value) {
            data.best_value = node.value + lb;
        }
        if (node.value + ub > data.best_value) {
            /* pq.push({node.value + ub, node.index + 1, hint, node.free_capacity, node.value}); */
            pq.push({node.value + ub, node.index + 1, node.free_capacity, node.value});
        }
        if (node.free_capacity >= instance.items()[node.index].weight) {
            node.value += instance.items()[node.index].value;
            node.free_capacity -= instance.items()[node.index].weight;
            ++node.index;
            pq.push(node);
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
    /* Node node{0, 0, 1, instance.capacity(), 0}; */
    Node node{0, 0, instance.capacity(), 0};
    std::clog << "Working...\n";
    auto t_start = std::chrono::steady_clock::now();
    /* auto const& [lb, ub] = instance.compute_bounds_hint(node.free_capacity, node.index, node.hint); */
    auto const& [lb, ub] = instance.compute_bounds_linear(node.free_capacity, node.index);
    data.best_value = lb;
    if (lb < ub) {
        node.upper_bound = ub;
        pq.push(node);
        knapsack(pq, data, instance);
    }
    auto t_end = std::chrono::steady_clock::now();
    std::clog << "Finished\n" << std::endl;

    std::clog << "Time (s): " << std::fixed << std::setprecision(3) << std::chrono::duration<double>(t_end - t_start).count() << '\n';
    std::clog << "Processed nodes: " << data.processed_nodes << '\n';
    std::clog << "Solution: " << data.best_value << '\n';

    std::cout << "time,processed_nodes,solution\n";
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count() << ','
              << data.processed_nodes << ',' << data.best_value << std::endl;
}
