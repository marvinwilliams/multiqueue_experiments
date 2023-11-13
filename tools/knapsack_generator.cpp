#include "cxxopts.hpp"

#include <iomanip>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

template <typename WeightType = long long, typename ValueType = WeightType>
struct Item {
    WeightType weight;
    ValueType value;
};

template <typename WeightType = long long, typename ValueType = WeightType>
void generate_knapsack_instance(long long n, WeightType a, WeightType b, ValueType l, ValueType u, double f,
                                unsigned long s) {
    std::vector<Item<WeightType, ValueType>> elements(static_cast<std::size_t>(n));
    std::default_random_engine rng(s);
    std::generate(elements.begin(), elements.end(), [&rng, a, b, l, u]() {
        if constexpr (std::is_floating_point_v<WeightType>) {
            std::uniform_real_distribution<WeightType> random_weight(a, b);
            std::uniform_real_distribution<ValueType> random_add(l, u);
            auto weight = random_weight(rng);
            return Item<WeightType, ValueType>{weight, static_cast<ValueType>(weight) + random_add(rng)};
        } else {
            std::uniform_int_distribution<WeightType> random_weight(a, b);
            std::uniform_int_distribution<ValueType> random_add(l, u);
            auto weight = random_weight(rng);
            return Item<WeightType, ValueType>{weight, static_cast<ValueType>(weight) + random_add(rng)};
        }
    });
    std::cout << std::fixed << std::setprecision(6);
    double c = n * (b - a) / f;
    std::cout << n << '\n' << c << "\n\n";
    for (auto e : elements) {
        std::cout << e.value << ' ' << e.weight << '\n';
    }
    std::cout << '\n';
}

int main(int argc, char* argv[]) {
    cxxopts::Options options("Knapsack generator", "Generate knapsack instances");
    long long n = 1000;
    double a = 1000;
    double b = 100000;
    double l = 10000;
    double u = 12500;
    double f = 2;
    unsigned long s = 1;
    bool use_doubles = false;
    // clang-format off
    options.add_options()
      ("n,num-elements", "Number of elements", cxxopts::value<long long>(n), "NUMBER")
      ("a,weight-min", "Min weight", cxxopts::value<double>(a), "NUMBER")
      ("b,weight-max", "Max weight", cxxopts::value<double>(b), "NUMBER")
      ("l,p-min", "Min add to profits", cxxopts::value<double>(l), "NUMBER")
      ("u,p-max", "Max add to profits", cxxopts::value<double>(u), "NUMBER")
      ("f,capacity-divide", "Fraction of weights", cxxopts::value<double>(f), "NUMBER")
      ("s,seed", "Seed", cxxopts::value<unsigned long>(s), "NUMBER")
      ("d,doubles", "Use doubles", cxxopts::value<bool>(use_doubles)->default_value("false"), "Use doubles")
      ("h,help", "Print this help");
    // clang-format on

    try {
        auto result = options.parse(argc, argv);
        if (result.count("help") > 0) {
            std::cerr << options.help() << std::endl;
            return 0;
        }
    } catch (cxxopts::OptionParseException const& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    if (use_doubles) {
        generate_knapsack_instance<double>(n, a, b, l, u, f, s);
    } else {
        generate_knapsack_instance<long long, long long>(n, a, b, l, u, f, s);
    }
    return 0;
}
