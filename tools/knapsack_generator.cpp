#include "cxxopts.hpp"

#include <iostream>
#include <random>
#include <utility>
#include <vector>

struct Element {
    std::size_t weight;
    std::size_t value;
};

int main(int argc, char* argv[]) {
    cxxopts::Options options("Knapsack generator",
                             "Generate knapsack instances");
    // clang-format off
    options.add_options()
      ("n,num-elements", "Number of elements", cxxopts::value<std::size_t>(), "NUMBER")
      ("a,weight-min", "Min weight", cxxopts::value<std::size_t>(), "NUMBER")
      ("b,weight-max", "Max weight", cxxopts::value<std::size_t>(), "NUMBER")
      ("l,p-min", "Min add to profits", cxxopts::value<std::size_t>(), "NUMBER")
      ("u,p-max", "Max add to profits", cxxopts::value<std::size_t>(), "NUMBER")
      ("f,capacity-divide", "Fraction of weights", cxxopts::value<std::size_t>(), "NUMBER")
      ("s,seed", "Seed", cxxopts::value<std::size_t>(), "NUMBER")
      ("h,help", "Print this help");
    // clang-format on

    std::size_t n = 1000;
    std::size_t a = 10;
    std::size_t b = 1000;
    std::size_t l = 0;
    std::size_t u = 200;
    std::size_t f = 2;
    std::size_t s = 1;

    try {
        auto result = options.parse(argc, argv);
        if (result.count("help")) {
            std::cerr << options.help() << std::endl;
            return 0;
        }
        if (result.count("num-elements") > 0) {
            n = result["num-elements"].as<std::size_t>();
        }
        if (result.count("weight-min") > 0) {
            a = result["weight-min"].as<std::size_t>();
        }
        if (result.count("weight-max") > 0) {
            b = result["weight-max"].as<std::size_t>();
        }
        if (result.count("p-min") > 0) {
            l = result["p-min"].as<std::size_t>();
        }
        if (result.count("p-max") > 0) {
            u = result["p-max"].as<std::size_t>();
        }
        if (result.count("capacity-divide") > 0) {
            f = result["capacity-divide"].as<std::size_t>();
        }
        if (result.count("seed") > 0) {
            s = result["seed"].as<std::size_t>();
        }
    } catch (cxxopts::OptionParseException const& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    std::vector<Element> elements(n);
    std::size_t sum = 0;
    std::mt19937_64 rng(s);
    std::uniform_int_distribution<std::size_t> random_weight(a, b);
    std::uniform_int_distribution<std::size_t> random_add(l, u);
    for (std::size_t i = 0; i < n; ++i) {
        elements[i].weight = random_weight(rng);
        elements[i].value = elements[i].weight + random_add(rng);
        sum += elements[i].weight;
    }
    std::size_t c = sum / f;
    std::cout << n << '\n' << c << "\n\n";
    for (auto e : elements) {
        std::cout << e.value << ' ' << e.weight << '\n';
    }
    std::cout << '\n';
}
