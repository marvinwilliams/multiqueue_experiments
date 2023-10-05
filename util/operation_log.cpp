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

// A push might be after the pop because the tick is recorded
// after/before the push/pop is performed.
bool operation_log::verify_logs(OperationLog const& logs) {
    auto popped = std::unordered_set<std::size_t>{};
    long long num_late_pushes = 0;
    for (auto op : logs.pops) {
        if (op.ref_index >= logs.pushes.size()) {
            std::cerr << "Error: Popped element has invalid index: " << op.ref_index << '\n';
            return false;
        }
        if (op.ref_index != logs.pushes[op.ref_index].index) {
            std::cerr << "Error: Pushes have wrong indices\n";
            return false;
        }
        if (op.tick < logs.pushes[op.ref_index].tick) {
            ++num_late_pushes;
        }
        if (!popped.insert(op.ref_index).second) {
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
    auto pop_it = logs.pops.begin();
    auto push_it = logs.pushes.begin();
    while (pop_it != logs.pops.end() && push_it != logs.pushes.end()) {
        if (pop_it->tick < push_it->tick) {
            out << "d," << pop_it->tick << ',' << pop_it->ref_index << '\n';
            ++pop_it;
        } else {
            out << "i," << push_it->tick << ',' << push_it->key << ',' << push_it->index << '\n';
            ++push_it;
        }
    }
    while (pop_it != logs.pops.end()) {
        out << "d," << pop_it->tick << ',' << pop_it->ref_index << '\n';
        ++pop_it;
    }
    while (push_it != logs.pushes.end()) {
        out << "i," << push_it->tick << ',' << push_it->key << ',' << push_it->index << '\n';
        ++push_it;
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

std::vector<operation_log::Metrics> operation_log::replay_logs(OperationLog const& logs) {
    auto push_lookup = std::vector<std::size_t>(logs.pushes.size());
    for (std::size_t i = 0; i < logs.pushes.size(); ++i) {
        push_lookup[logs.pushes[i].index] = i;
    }
    ReplayTree<unsigned long, HeapElement, ExtractKey> replay_tree{};
    std::vector<Metrics> metrics;
    metrics.reserve(logs.pops.size());
    auto push_it = logs.pushes.begin();
    for (auto const& pop : logs.pops) {
        auto const& push = logs.pushes[push_lookup[pop.ref_index]];
        // Inserting everything before next deletion
        auto until = std::max(push.tick, pop.tick);
        for (; push_it != logs.pushes.end() && push_it->tick <= until; ++push_it) {
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
