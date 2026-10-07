#pragma once

#include "replay_tree.hpp"

#include <algorithm>
#include <cstddef>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace quality {

struct Metrics {
    std::size_t rank_error;
    std::size_t delay;
};

struct Log {
    using key_type = long long;
    struct Pop {
        std::size_t pushes_before;
        std::size_t push_index;
    };
    std::vector<key_type> keys;
    std::vector<Pop> pops;

    [[nodiscard]] std::size_t invalid_pops() const noexcept {
        return static_cast<std::size_t>(std::count_if(
            pops.begin(), pops.end(), [](Pop const& p) { return p.push_index >= p.pushes_before; }));
    }
};

// Text format: "<pushes> <pops>", then one line per operation in order, "+<key>" or "-<push index>"
inline void write_log(Log const& log, std::ostream& out) {
    out << log.keys.size() << ' ' << log.pops.size() << '\n';
    std::size_t i = 0;
    for (auto const& pop : log.pops) {
        for (; i < pop.pushes_before; ++i) {
            out << '+' << log.keys[i] << '\n';
        }
        out << '-' << pop.push_index << '\n';
    }
    for (; i < log.keys.size(); ++i) {
        out << '+' << log.keys[i] << '\n';
    }
}

inline std::vector<Metrics> replay(Log const& log) {
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
    std::size_t inserted = 0;
    for (auto const& pop : log.pops) {
        if (pop.push_index >= log.keys.size()) {
            throw std::runtime_error{"Pop references nonexistent push " + std::to_string(pop.push_index)};
        }
        // A pop logged before its push is replayed as if the push happened right before it
        auto const insert_until = std::max(pop.pushes_before, pop.push_index + 1);
        for (; inserted < insert_until; ++inserted) {
            replay_tree.insert({log.keys[inserted], inserted});
        }
        auto [success, rank, delay] = replay_tree.erase_val({log.keys[pop.push_index], pop.push_index});
        if (!success) {
            throw std::runtime_error{"Failed to delete element " + std::to_string(pop.push_index) + " with key " +
                                     std::to_string(log.keys[pop.push_index])};
        }
        metrics.push_back({rank, delay});
    }
    return metrics;
}

struct Distribution {
    double mean = 0.0;
    std::size_t p50 = 0;
    std::size_t p90 = 0;
    std::size_t p99 = 0;
    std::size_t max = 0;

    template <typename JsonObject>
    void write_json(JsonObject& obj) const {
        obj.entry("mean", mean);
        obj.entry("p50", p50);
        obj.entry("p90", p90);
        obj.entry("p99", p99);
        obj.entry("max", max);
    }
};

inline Distribution distribution(std::vector<std::size_t> values) {
    Distribution d;
    if (values.empty()) {
        return d;
    }
    std::sort(values.begin(), values.end());
    auto quantile = [&values](double q) { return values[static_cast<std::size_t>(q * double(values.size() - 1))]; };
    double sum = 0.0;
    for (auto v : values) {
        sum += double(v);
    }
    d.mean = sum / double(values.size());
    d.p50 = quantile(0.5);
    d.p90 = quantile(0.9);
    d.p99 = quantile(0.99);
    d.max = values.back();
    return d;
}

struct Summary {
    Distribution rank_error;
    Distribution delay;

    template <typename JsonObject>
    void write_json(JsonObject& obj) const {
        obj.object("rank_error", [this](JsonObject& o) { rank_error.write_json(o); });
        obj.object("delay", [this](JsonObject& o) { delay.write_json(o); });
    }
};

inline Summary summarize(std::vector<Metrics> const& metrics) {
    std::vector<std::size_t> values(metrics.size());
    Summary s;
    std::transform(metrics.begin(), metrics.end(), values.begin(), [](Metrics const& m) { return m.rank_error; });
    s.rank_error = distribution(values);
    std::transform(metrics.begin(), metrics.end(), values.begin(), [](Metrics const& m) { return m.delay; });
    s.delay = distribution(std::move(values));
    return s;
}

}  // namespace quality
