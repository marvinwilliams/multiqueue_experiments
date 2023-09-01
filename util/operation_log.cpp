#include "operation_log.hpp"

#include "replay_tree.hpp"

#include <algorithm>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <istream>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <unordered_set>
#include <vector>

operation_log::OperationLog operation_log::merge_logs(std::vector<OperationLog> const& logs) {
    auto merged = OperationLog{};
    merged.pushes.resize(
        std::accumulate(logs.begin(), logs.end(), std::size_t{0},
                        [](std::size_t acc, OperationLog const& log) { return acc + log.pushes.size(); }));
    for (auto const& log : logs) {
        for (auto const& push : log.pushes) {
            merged.pushes[push.index] = push;
        }
    }
    merged.pops.reserve(
        std::accumulate(logs.begin(), logs.end(), std::size_t{0},
                        [](std::size_t acc, OperationLog const& log) { return acc + log.pops.size(); }));
    for (auto const& log : logs) {
        merged.pops.insert(merged.pops.end(), log.pops.begin(), log.pops.end());
    }
    std::sort(merged.pops.begin(), merged.pops.end());
    return merged;
}

// A push might be after the pop because the tick is recorded
// after/before the push/pop is performed.
bool operation_log::verify_logs(OperationLog const& logs) {
    auto popped = std::unordered_set<std::size_t>{};
    long long num_late_pushes = 0;
    for (auto op : logs.pops) {
        if (op.index >= logs.pushes.size()) {
            std::cerr << "Error: Popped element has invalid index: " << op.index << '\n';
            return false;
        }
        if (op.index != logs.pushes[op.index].index) {
            std::cerr << "Error: Pushes have wrong indices\n";
            return false;
        }
        if (op.tick < logs.pushes[op.index].tick) {
            ++num_late_pushes;
        }
        if (!popped.insert(op.index).second) {
            std::cerr << "Error: Popped the same element twice\n";
            return false;
        }
    }
    if (num_late_pushes > 0) {
        std::clog << "Warning: " << num_late_pushes << " elements inserted after their deletion\n";
    }
    return true;
}

void operation_log::write_logs(OperationLog const& logs, std::ostream& out) {
    out << logs.pushes.size() << ' ' << logs.pops.size() << '\n';
    for (auto const& push : logs.pushes) {
        out << "i," << push.tick << ',' << push.key << ',' << push.index << '\n';
    }
    for (auto const& pops : logs.pops) {
        out << "d," << pops.tick << ',' << pops.index << '\n';
    }
}

struct HeapElement {
    unsigned long key;
    std::size_t index;
    friend bool operator==(HeapElement const& lhs, HeapElement const& rhs) {
        return lhs.index == rhs.index;
    }
    friend bool operator!=(HeapElement const& lhs, HeapElement const& rhs) {
        return !(lhs == rhs);
    }
};

struct ExtractKey {
    static auto const& get(HeapElement const& e) {
        return e.key;
    }
};

std::vector<operation_log::Metrics> operation_log::replay_logs(OperationLog logs) {
    auto sorted_pushes = logs.pushes;
    std::sort(sorted_pushes.begin(), sorted_pushes.end());
    ReplayTree<unsigned long, HeapElement, ExtractKey> replay_tree{};
    auto push_it = sorted_pushes.begin();
    std::vector<Metrics> metrics;
    metrics.reserve(logs.pops.size());
    for (auto const& pop : logs.pops) {
        auto const& push = logs.pushes[pop.index];
        // Inserting everything before next deletion
        auto until = std::max(push.tick, pop.tick);
        for (; push_it != sorted_pushes.end() && push_it->tick <= until; ++push_it) {
            replay_tree.insert({push_it->key, push_it->index});
        }
        auto [success, rank, delay] = replay_tree.erase_val({push.key, push.index});
        if (!success) {
            std::cerr << "Failed to delete element " << push.index << " with key " << push.key << '\n';
            std::abort();
        }
        metrics.push_back({rank, delay});
    }
    return metrics;
}

void operation_log::write_metrics(std::vector<Metrics> const& metrics, std::ostream& out) {
    assert(!metrics.empty());
    for (auto const& m : metrics) {
        out << m.rank_error << ' ' << m.delay << '\n';
    }
}

void operation_log::write_metrics_average(std::vector<Metrics> const& metrics, std::ostream& out) {
    assert(!metrics.empty());
    auto summed_metrics = std::reduce(metrics.begin(), metrics.end(), Metrics{}, [](Metrics a, Metrics const& b) {
        a.rank_error += b.rank_error;
        a.delay += b.delay;
        return a;
    });
    out << "rank_error,delay\n";
    out << std::fixed << std::setprecision(2)
        << static_cast<double>(summed_metrics.rank_error) / static_cast<double>(metrics.size()) << ','
        << static_cast<double>(summed_metrics.delay) / static_cast<double>(metrics.size()) << '\n';
}
