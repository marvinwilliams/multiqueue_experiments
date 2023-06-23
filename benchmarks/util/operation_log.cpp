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
#include <vector>

// Note that a push might be after the pop because the tick is recorded
// after the push is performed.
void operation_log::verify_logs(std::vector<operation_log::OperationLog> const& logs) {
    for (auto const& l : logs) {
        if (!std::is_sorted(l.pushes.begin(), l.pushes.end())) {
            throw std::runtime_error{"Inconsistent push order"};
        }
        if (!std::is_sorted(l.pops.begin(), l.pops.end())) {
            throw std::runtime_error{"Inconsistent pop order"};
        }
    }
    auto popped = std::vector<std::vector<bool>>();
    for (auto const& l : logs) {
        popped.emplace_back(l.pushes.size(), false);
    }
    for (auto const& l : logs) {
        for (auto op : l.pops) {
            auto thread_id = extract_thread_id(op.data);
            auto elem_id = extract_elem_id(op.data);
            if (thread_id < 0 || thread_id >= static_cast<int>(logs.size())) {
                throw std::runtime_error{"Popped element from invalid thread: " + std::to_string(thread_id)};
            }
            auto const& thread_pushes = logs[static_cast<std::size_t>(thread_id)].pushes;
            if (elem_id >= thread_pushes.size()) {
                throw std::runtime_error{"Popped element has invalid element id: " + std::to_string(elem_id)};
            }
            if (popped[static_cast<std::size_t>(thread_id)][elem_id]) {
                throw std::runtime_error{"Popped the same element twice"};
            }
            popped[static_cast<std::size_t>(thread_id)][elem_id] = true;
        }
    }
}

void operation_log::write_logs(std::vector<OperationLog> const& logs, std::ostream& out) {
    out << logs.size() << '\n';
    for (auto const& l : logs) {
        out << l.pushes.size() << ' ' << l.pops.size() << '\n';
    }
    for (std::size_t i = 0; i < logs.size(); ++i) {
        for (auto const& push : logs[i].pushes) {
            out << push.tick << " i " << i << ' ' << push.key << '\n';
        }
    }
    for (std::size_t i = 0; i < logs.size(); ++i) {
        for (auto const& pops : logs[i].pops) {
            out << pops.tick << " d " << i << ' ' << extract_thread_id(pops.data) << ' ' << extract_elem_id(pops.data)
                << '\n';
        }
    }
}

struct Element {
    unsigned long key;
    unsigned long data;
    friend bool operator==(Element const& lhs, Element const& rhs) {
        return lhs.data == rhs.data;
    }
    friend bool operator!=(Element const& lhs, Element const& rhs) {
        return !(lhs == rhs);
    }
};

operation_log::Histogram operation_log::to_histogram(std::vector<OperationLog> const& logs) {
    struct ExtractKey {
        static auto const& get(Element const& e) {
            return e.key;
        }
    };

    auto num_pops = std::accumulate(logs.begin(), logs.end(), 0UL,
                                    [](std::size_t sum, auto const& l) { return sum + l.pops.size(); });
    auto num_pushes = std::accumulate(logs.begin(), logs.end(), 0UL,
                                      [](std::size_t sum, auto const& l) { return sum + l.pushes.size(); });

    Histogram histogram{};
    histogram.ranks.reserve(num_pushes);
    histogram.delays.reserve(num_pops);

    ReplayTree<unsigned long, Element, ExtractKey> replay_tree{};
    std::vector<std::size_t> push_indices(logs.size(), 0);
    std::vector<std::size_t> pop_indices(logs.size(), 0);
    auto next_pop = [&logs, &pop_indices]() {
        auto min_tick = std::numeric_limits<std::size_t>::max();
        auto min_thread = std::numeric_limits<std::size_t>::max();
        for (std::size_t i = 0; i < logs.size(); ++i) {
            if (pop_indices[i] < logs[i].pops.size() && logs[i].pops[pop_indices[i]].tick < min_tick) {
                min_tick = logs[i].pops[pop_indices[i]].tick;
                min_thread = i;
            }
        }
        return std::make_pair(min_tick, min_thread);
    };
    for (std::size_t i = 0; i < num_pops; ++i) {
        auto [tick, thread] = next_pop();
        assert(thread != std::numeric_limits<std::size_t>::max());
        auto const& pop = logs[thread].pops[pop_indices[thread]];
        auto push_thread = static_cast<std::size_t>(extract_thread_id(pop.data));
        auto elem_id = extract_elem_id(pop.data);
        auto key = logs[push_thread].pushes[elem_id].key;
        // Inserting everything before next deletion
        for (std::size_t j = 0; j < logs.size(); ++j) {
            auto const& p = logs[j].pushes;
            for (; push_indices[j] < logs[j].pushes.size() && p[push_indices[j]].tick <= tick; ++push_indices[j]) {
                replay_tree.insert({p[push_indices[j]].key, pack(static_cast<int>(j), push_indices[j])});
            }
        }
        for (; push_indices[push_thread] <= elem_id; ++push_indices[push_thread]) {
            replay_tree.insert({logs[push_thread].pushes[push_indices[push_thread]].key,
                                pack(static_cast<int>(push_thread), push_indices[push_thread])});
        }
        auto [success, rank, delay] = replay_tree.erase_val({key, pop.data});
        assert(success);
        if (histogram.ranks.size() <= rank) {
            histogram.ranks.resize(rank + 1);
        }
        ++histogram.ranks[rank];
        if (histogram.delays.size() <= delay) {
            histogram.delays.resize(delay + 1);
        }
        ++histogram.delays[delay];
        ++pop_indices[thread];
    }
    return histogram;
}

void operation_log::write_histogram(Histogram const& histogram, std::ostream& out) {
    std::size_t i = 0;
    for (; i < std::min(histogram.ranks.size(), histogram.delays.size()); ++i) {
        out << i << ' ' << histogram.ranks[i] << ' ' << histogram.delays[i] << '\n';
    }
    for (; i < histogram.ranks.size(); ++i) {
        out << i << ' ' << histogram.ranks[i] << ' ' << 0 << '\n';
    }
    for (; i < histogram.delays.size(); ++i) {
        out << i << ' ' << 0 << ' ' << histogram.delays[i] << '\n';
    }
}
