#include "evaluate_quality.hpp"
#include "replay_tree.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <istream>
#include <numeric>
#include <optional>
#include <vector>

struct HeapElement {
    key_type key;
    unsigned int from_thread_id;
    std::size_t op_id;
};

struct ExtractKey {
    static key_type const& get(HeapElement const& e) {
        return e.key;
    }
};

bool verify(PushLogType& push_log, PopLogType const& pop_log) {
    std::vector<std::vector<bool>> popped;
    for (auto const& i : push_log) {
        popped.emplace_back(i.size(), false);
    }
    for (auto const& thread_pops : pop_log) {
        for (auto const& entry : thread_pops) {
            if (!entry.payload) {
                continue;
            }
            if (entry.tick < push_log[entry.payload->from_thread_id][entry.payload->op_id].tick) {
                if (entry.payload->op_id > 0 &&
                    entry.tick < push_log[entry.payload->from_thread_id][entry.payload->op_id - 1].tick) {
                    std::cerr << "Popped element before push of previous element finished" << std::endl;
                    return false;
                }
                push_log[entry.payload->from_thread_id][entry.payload->op_id].tick = entry.tick;
            }

            if (popped[entry.payload->from_thread_id][entry.payload->op_id]) {
                std::cerr << "Popped element twice" << std::endl;
                return false;
            }
            popped[entry.payload->from_thread_id][entry.payload->op_id] = true;
        }
    }
    return true;
}

bool evaluate(PushLogType& push_log, PopLogType const& pop_log, std::filesystem::path const& rank_file,
              std::filesystem::path const& delay_file) {
    std::clog << "Verifying operations\n";
    if (!verify(push_log, pop_log)) {
        return false;
    }
    std::clog << "Sorting pops\n";
    std::vector<PopLogEntry> all_pops;
    all_pops.reserve(std::accumulate(pop_log.begin(), pop_log.end(), 0UL,
                                     [](std::size_t sum, auto const& p) { return sum + p.size(); }));
    for (auto const& p : pop_log) {
        all_pops.insert(all_pops.end(), p.begin(), p.end());
    }
    std::sort(all_pops.begin(), all_pops.end(), [](auto const& lhs, auto const& rhs) { return lhs.tick < rhs.tick; });

    std::vector<std::size_t> ranks;
    ranks.reserve(all_pops.size());
    std::vector<std::size_t> delays;
    delays.reserve(all_pops.size());

    std::clog << "Replaying operations";
    ReplayTree<unsigned long, HeapElement, ExtractKey> replay_tree{};
    std::vector<std::size_t> push_id(push_log.size());
    std::size_t progress_counter = 0;
    std::size_t progress_update = (all_pops.size() / 10UL);
    for (auto& pop : all_pops) {
        // Inserting everything before next deletion
        for (unsigned int t = 0; t < push_log.size(); ++t) {
            while (push_id[t] < push_log[t].size() && push_log[t][push_id[t]].tick <= pop.tick) {
                replay_tree.insert({push_log[t][push_id[t]].key, t, push_id[t]});
                ++push_id[t];
            }
        }

        if (!pop.payload) {
            ranks.push_back(replay_tree.size());
            replay_tree.increase_global_delay();
            continue;
        }

        // Check if the insertion corresponding to the current deletion is already inserted
        if (push_id[pop.payload->from_thread_id] < pop.payload->op_id) {
            std::cerr << "Push after pop after verifying, this must not happen" << std::endl;
            return false;
        }

        auto key = push_log[pop.payload->from_thread_id][pop.payload->op_id].key;
        auto [start, end] = replay_tree.equal_range(key);
        start = std::find_if(start, end, [&pop](auto const& e) {
            return std::tie(pop.payload->from_thread_id, pop.payload->op_id) == std::tie(e.from_thread_id, e.op_id);
        });
        if (start == end) {
            std::cerr << "Element not in the heap at deletion time" << std::endl;
            return false;
        }
        ranks.push_back(replay_tree.get_rank(key));
        auto [success, delay] = replay_tree.erase(start);
        if (!success || delay < 0) {
            std::cerr << "Inconsistent replay tree state, this must not happen" << std::endl;
            return false;
        }
        delays.push_back(static_cast<std::size_t>(delay));
        replay_tree.increase_delay(key);
        if (progress_counter == progress_update) {
            progress_counter = 0;
            std::clog << '.';
        }
        ++progress_counter;
    }

    std::clog << "\nWriting histograms" << std::endl;
    if (auto out = std::ofstream{rank_file}; out) {
        std::vector<std::size_t> rank_histogram(*std::max_element(ranks.begin(), ranks.end()) + 1);
        for (auto r : ranks) {
            ++rank_histogram[r];
        }
        for (size_t i = 0; i < rank_histogram.size(); ++i) {
            out << i << " " << rank_histogram[i] << '\n';
        }
    } else {
        std::cerr << "Failed to open file to write rank histogram\n";
    }
    if (auto out = std::ofstream{delay_file}; out) {
        std::vector<std::size_t> delay_histogram(*std::max_element(delays.begin(), delays.end()) + 1);
        for (auto r : delays) {
            ++delay_histogram[r];
        }
        for (size_t i = 0; i < delay_histogram.size(); ++i) {
            out << i << " " << delay_histogram[i] << '\n';
        }
    } else {
        std::cerr << "Failed to open file to write delay histogram\n";
    }
    std::clog << "Done" << std::endl;
    return true;
}
