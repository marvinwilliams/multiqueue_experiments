#include "replay_tree.hpp"

#include <cstddef>
#include <iostream>
#include <vector>

struct Metrics {
    std::size_t rank_error;
    std::size_t delay;
};

void write_metrics(std::ostream& out, std::vector<Metrics> const& metrics) {
    out << "rank_error,delay\n";
    for (auto const& m : metrics) {
        out << m.rank_error << ',' << m.delay << '\n';
    }
}

struct Log {
    using key_type = long long;
    struct Pop {
        std::size_t push_index;
        std::size_t ref_index;
    };
    std::vector<key_type> keys;
    std::vector<Pop> pops;
};

Log read_log(std::istream& in) {
    std::size_t num_pushes = 0;
    std::size_t num_pops = 0;
    long long invalid_pops = 0;
    in >> num_pushes >> num_pops;
    Log log;
    log.keys.reserve(num_pushes);
    log.pops.reserve(num_pops);
    for (std::size_t i = 0; i < num_pops; ++i) {
        in.ignore();
        while (in.get() == '+') {
            Log::key_type key;  // NOLINT
            in >> key;
            log.keys.push_back(key);
            in.ignore();
        }
        std::size_t index;  // NOLINT
        in >> index;
        std::size_t push_index = log.keys.size();
        if (index >= log.keys.size()) {
            ++invalid_pops;
            push_index = index + 1;
        }
        log.pops.push_back({push_index, index});
    }
    std::cerr << "Invalid pops: " << invalid_pops << '\n';
    return log;
}

std::vector<Metrics> replay(Log const& log) {
    struct HeapElement {
        Log::key_type key;
        std::size_t index;
        bool operator==(HeapElement const& other) const noexcept {
            return index == other.index;
        }
        bool operator!=(HeapElement const& other) const noexcept {
            return !(*this == other);
        }
    };

    struct ExtractKey {
        static auto const& get(HeapElement const& e) noexcept {
            return e.key;
        }
    };

    ReplayTree<Log::key_type, HeapElement, ExtractKey> replay_tree{};
    std::vector<Metrics> metrics;
    metrics.reserve(log.pops.size());
    std::size_t push_index = 0;
    for (auto const& pop : log.pops) {
        for (; push_index < pop.push_index; ++push_index) {
            replay_tree.insert({log.keys[push_index], push_index});
        }
        auto [success, rank, delay] = replay_tree.erase_val({log.keys[pop.ref_index], pop.ref_index});
        if (!success) {
            std::cerr << "Failed to delete element " << pop.ref_index << " with key " << log.keys[pop.ref_index]
                      << '\n';
            std::abort();
        }
        metrics.push_back({rank, delay});
    }
    return metrics;
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);
    std::clog << "Reading log...\n";
    auto log = read_log(std::cin);
    std::clog << "Analyzing log...\n";
    auto metrics = replay(log);
    std::clog << "Writing metrics...\n";
    write_metrics(std::cout, metrics);
}
