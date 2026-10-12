#pragma once

#include "replay_tree.hpp"

#include <algorithm>
#include <cstddef>
#include <iterator>
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
        std::size_t thread = 0;
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

struct Bin {
    std::size_t lo;
    std::size_t hi;
    std::size_t count;
};

// Exact bins below 8, then 4 bins per power of two
inline std::size_t bin_shift(std::size_t value) noexcept {
    std::size_t shift = 0;
    while ((value >> shift) >= 8) {
        ++shift;
    }
    return shift;
}

inline std::vector<Bin> histogram(std::vector<std::size_t> values) {
    std::sort(values.begin(), values.end());
    std::vector<Bin> bins;
    for (auto v : values) {
        if (bins.empty() || v >= bins.back().hi) {
            auto shift = bin_shift(v);
            auto lo = (v >> shift) << shift;
            bins.push_back({lo, lo + (std::size_t{1} << shift), 0});
        }
        ++bins.back().count;
    }
    return bins;
}

struct Window {
    std::size_t first_pop;
    std::size_t pops;
    Distribution rank_error;
    Distribution delay;
};

struct Trace {
    std::size_t thread = 0;
    std::size_t first = 0;
    std::vector<Metrics> metrics;
};

struct Summary {
    Distribution rank_error;
    Distribution delay;
    std::vector<Bin> rank_error_histogram;
    std::vector<Bin> delay_histogram;
    std::vector<Window> history;
    Trace trace;

    template <typename JsonObject>
    void write_json(JsonObject& obj) const {
        auto write_bin = [](std::ostream& out, Bin const& b) { out << '[' << b.lo << ',' << b.hi << ',' << b.count << ']'; };
        obj.object("rank_error", [&](JsonObject& o) {
            rank_error.write_json(o);
            o.array("histogram", rank_error_histogram.begin(), rank_error_histogram.end(), write_bin);
        });
        obj.object("delay", [&](JsonObject& o) {
            delay.write_json(o);
            o.array("histogram", delay_histogram.begin(), delay_histogram.end(), write_bin);
        });
        obj.array("history", history.begin(), history.end(), [](std::ostream& out, Window const& w) {
            JsonObject o{out};
            o.entry("first_pop", w.first_pop);
            o.entry("pops", w.pops);
            o.object("rank_error", [&w](JsonObject& d) { w.rank_error.write_json(d); });
            o.object("delay", [&w](JsonObject& d) { w.delay.write_json(d); });
        });
        obj.object("trace", [this](JsonObject& o) {
            o.entry("thread", trace.thread);
            o.entry("first", trace.first);
            o.array("rank_error", trace.metrics.begin(), trace.metrics.end(),
                    [](std::ostream& out, Metrics const& m) { out << m.rank_error; });
            o.array("delay", trace.metrics.begin(), trace.metrics.end(),
                    [](std::ostream& out, Metrics const& m) { out << m.delay; });
        });
    }
};

inline std::vector<std::size_t> rank_errors(std::vector<Metrics>::const_iterator first,
                                            std::vector<Metrics>::const_iterator last) {
    std::vector<std::size_t> values;
    values.reserve(static_cast<std::size_t>(last - first));
    std::transform(first, last, std::back_inserter(values), [](Metrics const& m) { return m.rank_error; });
    return values;
}

inline std::vector<std::size_t> delays(std::vector<Metrics>::const_iterator first,
                                       std::vector<Metrics>::const_iterator last) {
    std::vector<std::size_t> values;
    values.reserve(static_cast<std::size_t>(last - first));
    std::transform(first, last, std::back_inserter(values), [](Metrics const& m) { return m.delay; });
    return values;
}

inline std::vector<Window> history(std::vector<Metrics> const& metrics, std::size_t num_windows) {
    std::vector<Window> windows;
    num_windows = std::min(num_windows, metrics.size());
    for (std::size_t w = 0; w < num_windows; ++w) {
        auto first = metrics.size() * w / num_windows;
        auto last = metrics.size() * (w + 1) / num_windows;
        auto begin = metrics.begin() + static_cast<std::ptrdiff_t>(first);
        auto end = metrics.begin() + static_cast<std::ptrdiff_t>(last);
        windows.push_back({first, last - first, distribution(rank_errors(begin, end)), distribution(delays(begin, end))});
    }
    return windows;
}

inline Trace trace(Log const& log, std::vector<Metrics> const& metrics, std::size_t thread, std::size_t length) {
    Trace t;
    t.thread = thread;
    std::vector<std::size_t> indices;
    for (std::size_t i = 0; i < log.pops.size(); ++i) {
        if (log.pops[i].thread == thread) {
            indices.push_back(i);
        }
    }
    length = std::min(length, indices.size());
    t.first = (indices.size() - length) / 2;
    for (std::size_t i = t.first; i < t.first + length; ++i) {
        t.metrics.push_back(metrics[indices[i]]);
    }
    return t;
}

inline Summary summarize(std::vector<Metrics> const& metrics, std::size_t num_windows = 0) {
    Summary s;
    auto values = rank_errors(metrics.begin(), metrics.end());
    s.rank_error = distribution(values);
    s.rank_error_histogram = histogram(std::move(values));
    values = delays(metrics.begin(), metrics.end());
    s.delay = distribution(values);
    s.delay_histogram = histogram(std::move(values));
    s.history = history(metrics, num_windows);
    return s;
}

}  // namespace quality
