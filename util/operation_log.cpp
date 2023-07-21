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

// Note that a push might be after the pop because the tick is recorded
// after the push is performed.
bool operation_log::verify_logs(OperationLog const& logs) {
    auto popped = std::unordered_set<std::size_t>{};
    long long num_late_pushes = 0;
    for (auto op : logs.pops) {
        if (op.index >= logs.pushes.size()) {
            std::cerr << "Popped element has invalid index: " << op.index << '\n';
            return false;
        }
        if (op.tick < logs.pushes[op.index].tick) {
            ++num_late_pushes;
        }
        if (!popped.insert(op.index).second) {
            std::cerr << "Popped the same element twice\n";
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
        out << "i," << push.tick << ',' << push.key << '\n';
    }
    for (auto const& pops : logs.pops) {
        out << "d," << pops.tick << ',' << pops.index << '\n';
    }
}

struct Element {
    unsigned long key;
    std::size_t index;
    friend bool operator==(Element const& lhs, Element const& rhs) {
        return lhs.index == rhs.index;
    }
    friend bool operator!=(Element const& lhs, Element const& rhs) {
        return !(lhs == rhs);
    }
};

std::vector<operation_log::Metrics> operation_log::replay_logs(OperationLog logs) {
    struct ExtractKey {
        static auto const& get(Element const& e) {
            return e.key;
        }
    };
    std::vector<std::pair<Push, std::size_t>> sorted_pushes(logs.pushes.size());
    for (std::size_t i = 0; i < logs.pushes.size(); ++i) {
        sorted_pushes[i] = {logs.pushes[i], i};
    }
    std::sort(sorted_pushes.begin(), sorted_pushes.end(),
              [](auto const& lhs, auto const& rhs) { return lhs.first < rhs.first; });
    std::sort(logs.pops.begin(), logs.pops.end());

    ReplayTree<unsigned long, Element, ExtractKey> replay_tree{};
    auto push_it = sorted_pushes.begin();
    std::vector<Metrics> metrics;
    metrics.reserve(logs.pops.size());
    for (auto const& p : logs.pops) {
        auto key = logs.pushes[p.index].key;
        // Inserting everything before next deletion
        auto until = std::max(logs.pushes[p.index].tick, p.tick);
        for (; push_it != sorted_pushes.end() && push_it->first.tick <= until; ++push_it) {
            replay_tree.insert({push_it->first.key, push_it->second});
        }
        auto [success, rank, delay] = replay_tree.erase_val({key, p.index});
        if (!success) {
            std::cerr << "Failed to delete element " << p.index << " with key " << key << '\n';
            std::abort();
        }
        metrics.push_back({rank, delay});
    }
    return metrics;
}

void operation_log::write_metrics(std::vector<Metrics> const& metrics, std::ostream& out) {
    for (auto const& m : metrics) {
        out << m.rank_error << ' ' << m.delay << '\n';
    }
}
