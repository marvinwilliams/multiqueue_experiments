#include "cxxopts.hpp"

#include <iostream>
#include <random>
#include <utility>
#include <vector>

struct Element {
    long long weight;
    long long value;
};

int main(int argc, char* argv[]) {
    cxxopts::Options options("Knapsack generator", "Generate knapsack instances");
    // clang-format off
    options.add_options()
      ("n,num-elements", "Number of elements", cxxopts::value<long long>(), "NUMBER")
      ("a,weight-min", "Min weight", cxxopts::value<long long>(), "NUMBER")
      ("b,weight-max", "Max weight", cxxopts::value<long long>(), "NUMBER")
      ("l,p-min", "Min add to profits", cxxopts::value<long long>(), "NUMBER")
      ("u,p-max", "Max add to profits", cxxopts::value<long long>(), "NUMBER")
      ("f,capacity-divide", "Fraction of weights", cxxopts::value<long long>(), "NUMBER")
      ("s,seed", "Seed", cxxopts::value<std::size_t>(), "NUMBER")
      ("h,help", "Print this help");
    // clang-format on

    long long n = 1000;
    long long a = 1000;
    long long b = 100000;
    long long l = 10000;
    long long u = 12500;
    long long f = 2;
    std::size_t s = 1;

    try {
        auto result = options.parse(argc, argv);
        if (result.count("help") > 0) {
            std::cerr << options.help() << std::endl;
            return 0;
        }
        if (result.count("num-elements") > 0) {
            n = result["num-elements"].as<long long>();
        }
        if (result.count("weight-min") > 0) {
            a = result["weight-min"].as<long long>();
        }
        if (result.count("weight-max") > 0) {
            b = result["weight-max"].as<long long>();
        }
        if (result.count("p-min") > 0) {
            l = result["p-min"].as<long long>();
        }
        if (result.count("p-max") > 0) {
            u = result["p-max"].as<long long>();
        }
        if (result.count("capacity-divide") > 0) {
            f = result["capacity-divide"].as<long long>();
        }
        if (result.count("seed") > 0) {
            s = result["seed"].as<std::size_t>();
        }
    } catch (cxxopts::OptionParseException const& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    std::vector<Element> elements(static_cast<std::size_t>(n));
    std::default_random_engine rng(s);
    std::generate(elements.begin(), elements.end(), [&rng, a, b, l, u]() {
        std::uniform_int_distribution<long long> random_weight(a, b);
        std::uniform_int_distribution<long long> random_add(l, u);
        auto weight = random_weight(rng);
        return Element{weight, weight + random_add(rng)};
    });
    long long c = n * (b - a) / f;
    std::cout << n << '\n' << c << "\n\n";
    for (auto e : elements) {
        std::cout << e.value << ' ' << e.weight << '\n';
    }
    std::cout << '\n';
}
