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

operation_log::SortedLogs operation_log::sort_logs(std::vector<OperationLog> const& logs) {
    auto num_pops = std::accumulate(logs.begin(), logs.end(), 0UL,
                                    [](std::size_t sum, auto const& l) { return sum + l.pops.size(); });
    auto num_pushes = std::accumulate(logs.begin(), logs.end(), 0UL,
                                      [](std::size_t sum, auto const& l) { return sum + l.pushes.size(); });
    SortedLogs sorted_logs;
    sorted_logs.pushes.reserve(num_pushes);
    sorted_logs.pops.reserve(num_pops);

    std::vector<std::vector<bool>> popped;
    popped.reserve(logs.size());
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
            if (thread_pushes[elem_id].key != op.key) {
                throw std::runtime_error{"Popped element has wrong key"};
            }
            if (op.tick < thread_pushes[elem_id].tick) {
                op.tick = thread_pushes[elem_id].tick;
            }

            if (popped[static_cast<std::size_t>(thread_id)][elem_id]) {
                throw std::runtime_error{"Popped the same element twice"};
            }
            popped[static_cast<std::size_t>(thread_id)][elem_id] = true;
            sorted_logs.pops.push_back(op);
        }
    }
    std::sort(sorted_logs.pops.begin(), sorted_logs.pops.end(),
              [](auto const& lhs, auto const& rhs) { return lhs.tick < rhs.tick; });
    for (int i = 0; i < static_cast<int>(logs.size()); ++i) {
        for (std::size_t id = 0; id < logs[static_cast<std::size_t>(i)].pushes.size(); ++id) {
            auto const& e = logs[static_cast<std::size_t>(i)].pushes[id];
            sorted_logs.pushes.push_back({e.tick, {e.key, pack(i, id)}});
        }
    }
    std::sort(sorted_logs.pushes.begin(), sorted_logs.pushes.end(),
              [](auto const& lhs, auto const& rhs) { return lhs.tick < rhs.tick; });
    return sorted_logs;
}

operation_log::Histogram operation_log::to_histogram(SortedLogs const& logs) {
    struct ExtractKey {
        static unsigned long get(Element const& e) {
            return e.key;
        }
    };

    std::vector<std::size_t> ranks;
    ranks.reserve(logs.pops.size());
    std::size_t max_rank = 0;
    std::vector<std::size_t> delays;
    delays.reserve(logs.pops.size());
    std::size_t max_delay = 0;

    ReplayTree<unsigned long, Element, ExtractKey> replay_tree{};
    auto it = logs.pushes.begin();
    for (auto const& pop : logs.pops) {
        // Inserting everything before next deletion
        while (it != logs.pushes.end() && it->tick < pop.tick) {
            replay_tree.insert({it->element.key, it->element.data});
            ++it;
        }
        auto [start, end] = replay_tree.equal_range(pop.key);
        start = std::find_if(start, end, [&pop](auto const& e) { return e.data == pop.data; });
        assert(start != end);
        ranks.push_back(replay_tree.get_rank(pop.key));
        if (ranks.back() > max_rank) {
            max_rank = ranks.back();
        }
        auto [success, delay] = replay_tree.erase(start);
        assert(success && delay >= 0);
        delays.push_back(static_cast<std::size_t>(delay));
        if (delays.back() > max_delay) {
            max_delay = delays.back();
        }
        replay_tree.increase_delay(pop.key);
    }
    Histogram histogram;
    histogram.ranks.resize(max_rank + 1);
    histogram.delays.resize(max_delay + 1);
    for (auto r : ranks) {
        ++histogram.ranks[r];
    }
    for (auto d : delays) {
        ++histogram.delays[d];
    }
    return histogram;
}

void operation_log::write_logs(SortedLogs const& logs, std::filesystem::path const& file) {
    std::ofstream out(file);
    if (!out) {
        throw std::runtime_error("Failed to open file: " + file.string());
    }
    auto it = logs.pushes.begin();
    for (auto const& pop : logs.pops) {
        while (it != logs.pushes.end() && it->tick < pop.tick) {
            out << it->tick << "i " << extract_thread_id(it->element.data) << ' ' << it->element.key << ' '
                << it->element.data << '\n';
        }
        out << pop.tick << "d " << extract_thread_id(it->element.data) << pop.key << ' ' << it->element.data << '\n';
    }
    while (it != logs.pushes.end()) {
        out << it->tick << "i " << extract_thread_id(it->element.data) << ' ' << it->element.key << ' '
            << it->element.data << '\n';
    }
}
