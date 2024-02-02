#include <cxxopts.hpp>

#include <iomanip>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

void generate_uncorrelated(int n, int min_weight, int max_weight, int min_value, int max_value, double factor,
                           unsigned long seed) {
    std::cout << n << '\n';
    std::cout << static_cast<int>(static_cast<double>(n) *
                                  (static_cast<double>(min_weight) + static_cast<double>(max_weight - min_weight) / 2) *
                                  factor)
              << '\n'
              << '\n';
    std::seed_seq seq{seed};
    std::mt19937 rng(seq);
    rng.discard(700'000);
    std::uniform_int_distribution<int> weight(min_weight, max_weight);
    std::uniform_int_distribution<int> value(min_value, max_value);
    for (int i = 0; i < n; ++i) {
        auto w = weight(rng);
        auto v = value(rng);
        std::cout << v << " " << w << '\n';
    }
}

void generate_uncorrelated(int n, double min_weight, double max_weight, double min_value, double max_value,
                           double factor, unsigned long seed) {
    std::cout << n << '\n';
    std::cout << static_cast<int>(static_cast<double>(n) * (min_weight + (max_weight - min_weight) / 2) * factor)
              << '\n'
              << '\n';
    std::seed_seq seq{seed};
    std::mt19937 rng(seq);
    rng.discard(700'000);
    std::uniform_real_distribution<double> weight(min_weight, max_weight);
    std::uniform_real_distribution<double> value(min_value, max_value);
    for (int i = 0; i < n; ++i) {
        auto w = weight(rng);
        auto v = value(rng);
        std::cout << v << " " << w << '\n';
    }
}

void generate_correlated(int n, int min_weight, int max_weight, int min_value, int max_value, double factor,
                         unsigned long seed) {
    std::cout << n << '\n';
    std::cout << static_cast<int>(static_cast<double>(n) *
                                  (static_cast<double>(min_weight) + static_cast<double>(max_weight - min_weight) / 2) *
                                  factor)
              << '\n'
              << '\n';
    std::seed_seq seq{seed};
    std::mt19937 rng(seq);
    rng.discard(700'000);
    std::uniform_int_distribution<int> weight(min_weight, max_weight);
    std::uniform_int_distribution<int> value(min_value, max_value);
    for (int i = 0; i < n; ++i) {
        auto w = weight(rng);
        int v{};
        do {
            v = w + value(rng);
        } while (v < 0);
        std::cout << v << " " << w << '\n';
    }
}

void generate_correlated(int n, double min_weight, double max_weight, double min_value, double max_value, double factor,
                         unsigned long seed) {
    std::cout << n << '\n';
    std::cout << static_cast<int>(static_cast<double>(n) * (min_weight + (max_weight - min_weight) / 2) * factor)
              << '\n'
              << '\n';
    std::seed_seq seq{seed};
    std::mt19937 rng(seq);
    rng.discard(700'000);
    std::uniform_real_distribution<double> weight(min_weight, max_weight);
    std::uniform_real_distribution<double> value(min_value, max_value);
    for (int i = 0; i < n; ++i) {
        auto w = weight(rng);
        double v{};
        do {
            v = w + value(rng);
        } while (v < 0);
        std::cout << v << " " << w << '\n';
    }
}

int main(int argc, char* argv[]) {
    cxxopts::Options options("Knapsack generator", "Generate knapsack instances");
    int type = 0;
    int n = 1000;
    double min_weight = 1;
    double max_weight = 1000;
    double min_value = 0;
    double max_value = 100;
    double factor = 0.5;
    unsigned long seed = 1;
    // clang-format off
    options.add_options()
      ("t,type", "Type of instance (0: uncorrelated int, 1: correlated int, 2: uncorrelated real, 3: correlated real)", cxxopts::value<int>(type), "NUMBER")
      ("n,num-elements", "Number of elements", cxxopts::value<int>(n), "NUMBER")
      ("w,min-weight", "Minimum weight", cxxopts::value<double>(min_weight), "NUMBER")
      ("W,max-weight", "Maximum weight", cxxopts::value<double>(max_weight), "NUMBER")
      ("v,min-value", "Minimum value to profits", cxxopts::value<double>(min_value), "NUMBER")
      ("V,max-value", "Maximum value to profits", cxxopts::value<double>(max_value), "NUMBER")
      ("f,factor", "Capacity factor", cxxopts::value<double>(factor), "NUMBER")
      ("s,seed", "Seed", cxxopts::value<unsigned long>(seed), "NUMBER")
      ("h,help", "Print this help");
    // clang-format on

    try {
        auto result = options.parse(argc, argv);
        if (result.count("help") > 0) {
            std::cerr << options.help() << std::endl;
            return EXIT_SUCCESS;
        }
    } catch (cxxopts::OptionParseException const& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    if (type == 0) {
        generate_uncorrelated(n, static_cast<int>(min_weight), static_cast<int>(max_weight),
                              static_cast<int>(min_value), static_cast<int>(max_value), factor, seed);
    } else if (type == 1) {
        generate_correlated(n, static_cast<int>(min_weight), static_cast<int>(max_weight), static_cast<int>(min_value),
                            static_cast<int>(max_value), factor, seed);
    } else if (type == 2) {
        generate_uncorrelated(n, min_weight, max_weight, min_value, max_value, factor, seed);
    } else if (type == 3) {
        generate_correlated(n, min_weight, max_weight, min_value, max_value, factor, seed);
    } else {
        std::cerr << "Invalid type" << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
