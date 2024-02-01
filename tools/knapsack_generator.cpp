#include <cxxopts.hpp>

#include <iomanip>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

void generate_uncorrelated(int n, int max_weight, double factor, unsigned long seed) {
    std::cout << n << '\n';
    std::cout << static_cast<int>(static_cast<double>(n) * (static_cast<double>(max_weight) / 2 + 0.5) * factor)
              << '\n';
    std::default_random_engine rng(seed);
    std::uniform_int_distribution<int> weight(1, max_weight);
    std::uniform_int_distribution<int> value(0, max_weight);
    for (int i = 0; i < n; ++i) {
        std::cout << weight(rng) << " " << value(rng) << '\n';
    }
}

void generate_uncorrelated(int n, double max_weight, double factor, unsigned long seed) {
    std::cout << n << '\n';
    std::cout << static_cast<int>(static_cast<double>(n) * (static_cast<double>(max_weight) / 2 + 0.5) * factor)
              << '\n';
    std::default_random_engine rng(seed);
    std::uniform_real_distribution<double> weight(1.0, max_weight);
    std::uniform_real_distribution<double> value(0, max_weight);
    for (int i = 0; i < n; ++i) {
        std::cout << weight(rng) << " " << value(rng) << '\n';
    }
}

void generate_correlated(int n, int max_weight, int min_add, int max_add, double factor, unsigned long seed) {
    std::cout << n << '\n';
    std::cout << static_cast<int>(static_cast<double>(n) * (static_cast<double>(max_weight) / 2 + 0.5) * factor)
              << '\n';
    std::default_random_engine rng(seed);
    std::uniform_int_distribution<int> weight(1, max_weight);
    std::uniform_int_distribution<int> add(min_add, max_add);
    for (int i = 0; i < n; ++i) {
        auto w = weight(rng);
        std::cout << w << " " << w + add(rng) << '\n';
    }
}

void generate_correlated(int n, double max_weight, double min_add, double max_add, double factor,
                         unsigned long seed) {
    std::cout << n << '\n';
    std::cout << static_cast<int>(static_cast<double>(n) * (static_cast<double>(max_weight) / 2 + 0.5) * factor)
              << '\n';
    std::default_random_engine rng(seed);
    std::uniform_real_distribution<double> weight(1.0, max_weight);
    std::uniform_real_distribution<double> add(min_add, max_add);
    for (int i = 0; i < n; ++i) {
        auto w = weight(rng);
        std::cout << w << " " << w + add(rng) << '\n';
    }
}

int main(int argc, char* argv[]) {
    cxxopts::Options options("Knapsack generator", "Generate knapsack instances");
    int type = 0;
    int n = 1000;
    double max = 1000;
    double min_add = 100;
    double max_add = 500;
    double factor = 0.5;
    unsigned long seed = 1;
    // clang-format off
    options.add_options()
      ("t,type", "Type of instance (0: uncorrelated int, 1: correlated int, 2: uncorrelated real, 3: correlated real)", cxxopts::value<int>(type), "NUMBER")
      ("n,num-elements", "Number of elements", cxxopts::value<int>(n), "NUMBER")
      ("m,max-weight", "Maximum weight and value", cxxopts::value<double>(max), "NUMBER")
      ("l,min-add", "Minimum add to profits", cxxopts::value<double>(min_add), "NUMBER")
      ("u,max-add", "Maximum add to profits", cxxopts::value<double>(max_add), "NUMBER")
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
        generate_uncorrelated(n, static_cast<int>(max), factor, seed);
    } else if (type == 1) {
        generate_correlated(n, static_cast<int>(max), static_cast<int>(min_add), static_cast<int>(max_add), factor,
                            seed);
    } else if (type == 2) {
        generate_uncorrelated(n, max, factor, seed);
    } else if (type == 3) {
        generate_correlated(n, max, min_add, max_add, factor, seed);
    } else {
        std::cerr << "Invalid type" << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
