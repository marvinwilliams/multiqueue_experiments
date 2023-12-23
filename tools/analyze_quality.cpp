#include "benchmark/util/operation_log.hpp"
#include "util/replay_tree.hpp"

#include <cstddef>
#include <iostream>
#include <vector>

struct Metrics {
    std::size_t rank_error;
    std::size_t delay;
};

// Format:
// push_count,pop_count
// tick,key,index (repeat push_count times)
// tick,ref_index (repeat pop_count times)
std::OperationLog read(std::istream& in) {
  
}




    std::sort(op_log.pushes.begin(), op_log.pushes.end());
    std::sort(op_log.pops.begin(), op_log.pops.end());
    if (log_out.is_open()) {
        log(std::clog, context.start_time) << "Writing operation log...\n";
        operation_log::write(op_log, log_out);
        log_out.close();
    }
    log(std::clog, context.start_time) << "Replaying operations...\n";
    auto metrics = operation_log::replay(op_log);
    if (metrics_out.is_open()) {
        metrics_out << "rank_error,delay\n";
        for (auto const& m : metrics) {
            metrics_out << m.rank_error << ',' << m.delay << '\n';
        }
        metrics_out.close();
    }
    assert(!metrics.empty());
    auto summed_metrics =
        std::reduce(metrics.begin(), metrics.end(), operation_log::Metrics{}, [](auto a, auto const& b) {
            a.rank_error += b.rank_error;
            a.delay += b.delay;
            return a;
        });
std::vector<Metrics> replay(OperationLog const& logs) {
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

    struct HeapElement {
        unsigned long key;
        std::size_t index;
        bool operator==(HeapElement const& other) const noexcept {
            return index == other.index;
        }
        bool operator!=(HeapElement const& other) const noexcept {
            return !(*this == other);
        }
    };

    struct ExtractKey {
        static auto const& get(HeapElement const& e) {
            return e.key;
        }
    };

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
