#include "quality.hpp"

#include <cstddef>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <vector>

void write_metrics(std::ostream& out, std::vector<quality::Metrics> const& metrics) {
    out << "rank_error,delay\n";
    for (auto const& m : metrics) {
        out << m.rank_error << ',' << m.delay << '\n';
    }
}

quality::Log read_log(std::istream& in) {
    std::size_t num_pushes = 0;
    std::size_t num_pops = 0;
    in >> num_pushes >> num_pops;
    quality::Log log;
    log.keys.reserve(num_pushes);
    log.pops.reserve(num_pops);
    char op{};
    while (in >> op) {
        if (op == '+') {
            quality::Log::key_type key;  // NOLINT
            in >> key;
            log.keys.push_back(key);
        } else if (op == '-') {
            std::size_t index;  // NOLINT
            in >> index;
            log.pops.push_back({log.keys.size(), index});
        } else {
            std::cerr << "Invalid operation '" << op << "' in log\n";
            std::exit(EXIT_FAILURE);
        }
    }
    if (log.keys.size() != num_pushes || log.pops.size() != num_pops) {
        std::cerr << "Wrong number of pushes or pops\n";
        std::exit(EXIT_FAILURE);
    }
    std::cerr << "Invalid pops: " << log.invalid_pops() << '\n';
    return log;
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);
    std::clog << "Reading log...\n";
    auto log = read_log(std::cin);
    std::clog << "Analyzing log...\n";
    std::vector<quality::Metrics> metrics;
    try {
        metrics = quality::replay(log);
    } catch (std::exception const& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
    std::clog << "Writing metrics...\n";
    write_metrics(std::cout, metrics);
}
