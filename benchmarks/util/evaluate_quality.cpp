#include "evaluate_quality.hpp"

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

struct HeapElement {
    unsigned long key;
    int thread_id;
    std::size_t elem_id;
};

struct ExtractKey {
    static unsigned long const& get(HeapElement const& e) {
        return e.key;
    }
};

struct HistogramEntry {
    std::size_t rank_count;
    std::size_t delay_count;
};

using Histogram = std::vector<HistogramEntry>;

void quality::fix_logs(PushLogType& push_log, PopLogType const& pop_log) {
    std::vector<std::vector<bool>> popped;
    for (auto const& i : push_log) {
        popped.emplace_back(i.size(), false);
    }
    for (auto const& thread_pops : pop_log) {
        for (auto const& entry : thread_pops) {
            if (entry.thread_id < 0 || entry.thread_id >= static_cast<int>(push_log.size())) {
                throw std::runtime_error{"Popped element from invalid thread: " + std::to_string(entry.thread_id)};
            }
            auto& thread_push_log = push_log[static_cast<std::size_t>(entry.thread_id)];
            if (entry.elem_id >= thread_push_log.size()) {
                throw std::runtime_error{"Popped element has invalid element id: " + std::to_string(entry.elem_id)};
            }
            if (entry.tick < thread_push_log[entry.elem_id].tick) {
                if (entry.elem_id > 0 && entry.tick < thread_push_log[entry.elem_id - 1].tick) {
                    throw std::runtime_error{"Popped element before pushing thread pushed previous element"};
                }
                thread_push_log[entry.elem_id].tick = entry.tick;
            }

            if (popped[static_cast<std::size_t>(entry.thread_id)][entry.elem_id]) {
                throw std::runtime_error{"Popped element twice"};
            }
            popped[static_cast<std::size_t>(entry.thread_id)][entry.elem_id] = true;
        }
    }
}

void quality::write_histogram(PushLogType const& push_log, PopLogType const& pop_log,
                              std::filesystem::path const& file) {
    std::clog << "Preprocessing logs..." << std::flush;
    auto num_pops = std::accumulate(pop_log.begin(), pop_log.end(), 0UL,
                                    [](std::size_t sum, auto const& p) { return sum + p.size(); });
    std::vector<std::vector<PopLogEntry>::const_iterator> iters;
    iters.resize(pop_log.size());
    for (std::size_t t = 0; t < pop_log.size(); ++t) {
        iters[t] = pop_log[t].begin();
    }
    std::clog << "done" << std::endl;

    std::vector<std::size_t> ranks;
    ranks.reserve(num_pops);
    std::vector<std::size_t> delays;
    delays.reserve(num_pops);

    std::clog << "Replaying...";
    ReplayTree<unsigned long, HeapElement, ExtractKey> replay_tree{};
    std::vector<std::size_t> push_id(push_log.size());
    std::size_t max_entry = 0;
    for (std::size_t i = 0; i < num_pops; ++i) {
        auto min_it = iters.end();
        for (std::size_t i = 0; i < iters.size(); ++i) {
            if (iters[i] != pop_log[i].end()) {
                if (min_it == iters.end() || iters[i]->tick < (*min_it)->tick) {
                    min_it = iters.begin() + static_cast<std::ptrdiff_t>(i);
                }
            }
        }
        assert(min_it == iters.end());
        // Inserting everything before next deletion
        for (std::size_t t = 0; t < push_log.size(); ++t) {
            while (push_id[t] < push_log[t].size() && push_log[t][push_id[t]].tick <= (*min_it)->tick) {
                replay_tree.insert({push_log[t][push_id[t]].key, static_cast<int>(t), push_id[t]});
                ++push_id[t];
            }
        }

        // Check if the insertion corresponding to the current deletion is already inserted
        if (push_id[static_cast<std::size_t>((*min_it)->thread_id)] < (*min_it)->elem_id) {
            std::cerr << "Error: Element pushed after it was popped" << std::endl;
            std::abort();
        }

        auto key = push_log[static_cast<std::size_t>((*min_it)->thread_id)][(*min_it)->elem_id].key;
        auto [start, end] = replay_tree.equal_range(key);
        start = std::find_if(start, end, [&](auto const& e) {
            return std::tie((*min_it)->thread_id, (*min_it)->elem_id) == std::tie(e.thread_id, e.elem_id);
        });
        assert(start != end);
        ranks.push_back(replay_tree.get_rank(key));
        if (ranks.back() > max_entry) {
            max_entry = ranks.back();
        }
        auto [success, delay] = replay_tree.erase(start);
        assert(success && delay >= 0);
        delays.push_back(static_cast<std::size_t>(delay));
        if (delays.back() > max_entry) {
            max_entry = delays.back();
        }
        replay_tree.increase_delay(key);
        ++(*min_it);
    }
    Histogram histogram(max_entry + 1);
    for (auto r : ranks) {
        ++histogram[r].rank_count;
    }
    for (auto r : delays) {
        ++histogram[r].delay_count;
    }
    std::clog << "done\nWriting histogram..." << std::flush;
    std::ofstream out(file);
    if (!out) {
        throw std::runtime_error("Failed to open file: " + file.string());
    }
    out << num_pops << '\n';
    for (std::size_t i = 0; i < histogram.size(); ++i) {
        out << i << ' ' << histogram[i].rank_count << ' ' << histogram[i].delay_count << '\n';
    }
}

void quality::write_logs(PushLogType const& push_log, PopLogType const& pop_log, std::filesystem::path const& file) {
    std::ofstream out(file);
    if (!out) {
        throw std::runtime_error("Failed to open file: " + file.string());
    }
    std::vector<std::size_t> push_index(push_log.size());
    std::vector<std::size_t> pop_index(pop_log.size());
    enum class EntryType { Push, Pop, None };
    std::pair<EntryType, std::size_t> next_entry = {EntryType::None, 0};
    auto new_smallest = [&push_log, &pop_log, &push_index, &pop_index, &next_entry](auto key) {
        return next_entry.first == EntryType::None ||
            (key < (next_entry.first == EntryType::Push
                        ? push_log[next_entry.second][push_index[next_entry.second]].tick
                        : pop_log[next_entry.second][pop_index[next_entry.second]].tick));
    };
    while (true) {
        // find next entry with smallest tick
        next_entry = {EntryType::None, 0};
        for (std::size_t i = 0; i < push_index.size(); ++i) {
            if (push_index[i] < push_log[i].size() && new_smallest(push_log[i][push_index[i]].tick)) {
                next_entry = {EntryType::Push, i};
            }
        }
        for (std::size_t i = 0; i < pop_index.size(); ++i) {
            if (pop_index[i] < pop_log[i].size() && new_smallest(pop_log[i][pop_index[i]].tick)) {
                next_entry = {EntryType::Pop, i};
            }
        }
        // no more entries
        if (next_entry.first == EntryType::None) {
            break;
        }
        if (next_entry.first == EntryType::Push) {
            out << 'i' << ' ' << next_entry.second << ' '
                << push_log[next_entry.second][push_index[next_entry.second]].tick.time_since_epoch().count() << ' '
                << push_log[next_entry.second][push_index[next_entry.second]].key << '\n';
            ++push_index[next_entry.second];
        } else {
            out << 'd' << ' ' << next_entry.second << ' '
                << pop_log[next_entry.second][pop_index[next_entry.second]].tick.time_since_epoch().count() << ' '
                << pop_log[next_entry.second][pop_index[next_entry.second]].thread_id << ' '
                << pop_log[next_entry.second][pop_index[next_entry.second]].elem_id << '\n';
            ++pop_index[next_entry.second];
        }
    }
}
