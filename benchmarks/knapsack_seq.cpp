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

static constexpr auto default_n = 1000;
static constexpr auto default_f = 0.5;
#ifdef INTEGER_INSTANCE
using data_type = long long;
static constexpr data_type default_a = 1000;
static constexpr data_type default_b = 100000;
static constexpr data_type default_l = 10000;
static constexpr data_type default_u = 12500;
#else
using data_type = double;
static constexpr data_type default_a = 0.01;
static constexpr data_type default_b = 1.01;
static constexpr data_type default_l = 0.1;
static constexpr data_type default_u = 0.125;
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
    long long n = default_n;
    data_type a = default_a;
    data_type b = default_b;
    data_type l = default_l;
    data_type u = default_u;
    double f = default_f;
    unsigned long s = 1;

    void options(cxxopts::Options& cmd) {
        cmd.add_options()("n,num-elements", "Number of elements", cxxopts::value<long long>(n), "NUMBER")(
            "a,weight-min", "Min weight", cxxopts::value<data_type>(a), "NUMBER")(
            "b,weight-max", "Max weight", cxxopts::value<data_type>(b), "NUMBER")(
            "l,p-min", "Min add to profits", cxxopts::value<data_type>(l), "NUMBER")(
            "u,p-max", "Max add to profits", cxxopts::value<data_type>(u), "NUMBER")(
            "f,capacity-divide", "Fraction of weights", cxxopts::value<double>(f), "NUMBER")(
            "s,seed", "Seed", cxxopts::value<unsigned long>(s), "NUMBER");
    }
};

void write_settings(Settings const& settings, std::ostream& out) {
    out << "Number of items: " << settings.n << '\n'
        << "Instance parameters: "
        << "weight_min=" << settings.a << " weight_max=" << settings.b << " profit_add_min=" << settings.l
        << " profit_add_max=" << settings.u << " cap_fraction=" << settings.f << '\n'
        << "Seed: " << settings.s << '\n';
}

struct Data {
    data_type best_value{0};
    long long processed_nodes{0};
};

void knapsack(pq_type& pq, Data& data, KnapsackInstance<data_type> const& instance) noexcept {
    while (!pq.empty()) {
        auto node = pq.top();
        pq.pop();
        ++data.processed_nodes;
        if (node.upper_bound <= data.best_value) {
            return;
        }
        auto const& [lb, ub] = instance.compute_bounds_linear(node.free_capacity, node.index + 1);
        if (node.value + lb > data.best_value) {
            data.best_value = node.value + lb;
        }
        if (node.index + 2 < instance.size()) {
            if (node.value + ub > data.best_value) {
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
}

int main(int argc, char* argv[]) {
    write_build_info(std::clog);
    std::clog << "\nCommand line:";
    for (int i = 0; i < argc; ++i) {
        std::clog << ' ' << argv[i];
    }
    std::clog << '\n' << '\n';

    Settings settings{};
    cxxopts::Options cmd(argv[0]);
    cmd.add_options()("h,help", "Print this help");
    settings.options(cmd);

    try {
        auto result = cmd.parse(argc, argv);
        if (result.count("help") > 0U) {
            std::cerr << cmd.help() << std::endl;
            return EXIT_SUCCESS;
        }
    } catch (cxxopts::OptionParseException const& e) {
        std::cerr << "Error parsing arguments: " << e.what() << '\n';
        std::cerr << cmd.help() << std::endl;
        return EXIT_FAILURE;
    }

    write_settings(settings, std::clog);

    auto instance =
        KnapsackInstance<data_type>(settings.n, settings.a, settings.b, settings.l, settings.u, settings.f, settings.s);

    pq_type::container_type container;
    container.reserve(1 << 24);
    pq_type pq{{}, std::move(container)};
    Data data;
    Node node{0, 0, instance.capacity(), 0};
    std::clog << "Working...\n";
    auto t_start = std::chrono::steady_clock::now();
    auto const& [lb, ub] = instance.compute_bounds_linear(node.free_capacity, node.index);
    data.best_value = lb;
    if (lb < ub) {
        node.upper_bound = ub;
        pq.push(node);
        knapsack(pq, data, instance);
    }
    auto t_end = std::chrono::steady_clock::now();
    std::clog << "Finished\n" << std::endl;

    std::clog << "Time (s): " << std::fixed << std::setprecision(3)
              << std::chrono::duration<double>(t_end - t_start).count() << '\n';
    std::clog << "Processed nodes: " << data.processed_nodes << '\n';
    std::clog << "Solution: " << data.best_value << '\n';

    std::cout << "time,processed,solution\n";
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count() << ','
              << data.processed_nodes << ',' << data.best_value << std::endl;
}
