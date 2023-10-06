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

void operation_log::write(OperationLog const& logs, std::ostream& out) {
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

std::vector<operation_log::Metrics> operation_log::replay(OperationLog const& logs) {
    auto push_lookup = std::vector<std::size_t>(logs.pushes.size(), logs.pushes.size());
    for (std::size_t i = 0; i < logs.pushes.size(); ++i) {
        if (logs.pushes[i].index >= push_lookup.size()) {
            std::cerr << "Index " << logs.pushes[i].index << " is out of bounds\n";
            std::abort();
        }
        if (push_lookup[logs.pushes[i].index] != logs.pushes.size()) {
            std::cerr << "Index " << logs.pushes[i].index << " is not unique\n";
            std::abort();
        }
        push_lookup[logs.pushes[i].index] = i;
    }
    ReplayTree<unsigned long, HeapElement, ExtractKey> replay_tree{};
    std::vector<Metrics> metrics;
    metrics.reserve(logs.pops.size());
    auto push_it = logs.pushes.begin();
    long long wrong_order = 0;
    for (auto const& pop : logs.pops) {
        auto const& push = logs.pushes[push_lookup[pop.ref_index]];
        // Inserting everything before next deletion
        auto until_tick = pop.tick;
        if (push.tick > until_tick) {
            until_tick = push.tick;
            ++wrong_order;
        }
        for (; push_it != logs.pushes.end() && push_it->tick <= until_tick; ++push_it) {
            replay_tree.insert({push_it->key, push_it->index});
        }
        auto [success, rank, delay] = replay_tree.erase_val({push.key, push.index});
        if (!success) {
            std::cerr << "Failed to delete element " << push.index << " with key " << push.key << '\n';
            std::abort();
        }
        metrics.push_back({rank, delay});
    }
    if (wrong_order > 0) {
        std::cerr << "Warning: " << wrong_order << " elements were inserted after their deletion\n";
    }
    return metrics;
}
