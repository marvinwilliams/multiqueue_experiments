#pragma once

#include "util.hpp"

#include "cxxopts.hpp"

#include <cstddef>
#include <mutex>
#include <optional>
#include <queue>
#include <vector>

namespace wrapper::locked_pq {

template <bool Min = true, typename Key = unsigned long, typename Value = std::pair<unsigned long, unsigned long>,
          typename KeyOfValue = util::KeyOfValue<Key, Value>>
class LockedPQ {
   public:
    using key_type = Key;
    using value_type = Value;
    using key_compare = std::conditional_t<Min, std::greater<>, std::less<>>;
    using value_compare = util::ValueCompare<value_type, key_compare, KeyOfValue>;
    using handle_type = LockedPQ&;

   private:
    using pq_type = std::priority_queue<value_type, std::vector<value_type>, value_compare>;

    pq_type pq_{};
    std::mutex m_;

   public:
    LockedPQ(std::size_t initial_capacity) {
        std::vector<value_type> v{};
        v.reserve(initial_capacity);
        pq_ = pq_type{value_compare{}, std::move(v)};
    }

    void push(value_type const& value) {
        std::lock_guard<std::mutex> lock{m_};
        pq_.push(value);
    }

    std::optional<value_type> try_pop() {
        std::lock_guard<std::mutex> lock{m_};
        if (pq_.empty()) {
            return std::nullopt;
        }
        auto retval = pq_.top();
        pq_.pop();
        return retval;
    }

    handle_type get_handle() {
        return *this;
    }
};

template <bool Min = true, typename Key = unsigned long, typename Value = std::pair<unsigned long, unsigned long>,
          typename KeyOfValue = util::KeyOfValue<Key, Value>>
using PQWrapper = LockedPQ<Min, Key, Value, KeyOfValue>;

inline void add_options(cxxopts::Options& /*options*/) {
}

template <bool Min = true, typename Key = unsigned long, typename Value = std::pair<unsigned long, unsigned long>,
          typename KeyOfValue = util::KeyOfValue<Key, Value>>
LockedPQ<Min, Key, Value, KeyOfValue> create(int /*num_threads*/, std::size_t initial_capacity,
                                             cxxopts::ParseResult const& /*result*/) {
    return LockedPQ<Min, Key, Value, KeyOfValue>{initial_capacity};
}

template <typename PQ>
std::ostream& describe(PQ const& pq, std::ostream& out) {
    out << "Locked std::priority_queue";
    return out;
}

}  // namespace wrapper::locked_pq
