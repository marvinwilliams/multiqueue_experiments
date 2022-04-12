#include "cxxopts.hpp"

#include <math.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <thread>
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

static KnapsackInstance instance;

static value_type bestValue;

value_type getUpperBound(weight_type capacity, size_t index) {
    value_type result = 0;
    while (index < instance.items.size() &&
           capacity >= instance.items[index].weight) {
        result += instance.items[index].value;
        capacity -= instance.items[index].weight;
        ++index;
    }
    if (index < instance.items.size()) {
        double fraction =
            static_cast<double>(capacity) / instance.items[index].weight;
        result += std::ceil(fraction * instance.items[index].value);
    }
    return result;
}

void solveRec(weight_type capacity, size_t index, value_type currentValue) {
    if (currentValue > bestValue) {
        bestValue = currentValue;
    }

    if (capacity == 0 || index == instance.items.size()) {
        return;
    }

    value_type upperBound = getUpperBound(capacity, index);

    if (currentValue + upperBound <= bestValue) {
        return;
    }

    if (instance.items[index].weight > capacity) {
        solveRec(capacity, index + 1, currentValue);
        return;
    }

    solveRec(capacity - instance.items[index].weight, index + 1,
             currentValue + instance.items[index].value);
    solveRec(capacity, index + 1, currentValue);
    return;
}

static void read_problem(std::filesystem::path instance_file) {
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
    bestValue = 0;
    solveRec(instance.capacity, 0, 0);
    std::cout << "value: " << bestValue << '\n';
}
