#pragma once

#include "cxxopts.hpp"

#include <tbb/concurrent_priority_queue.h>

#include <optional>
#include <ostream>
#include <utility>

namespace wrapper {

template <bool Min, typename Key = unsigned long, typename T = unsigned long>
class TBBPriorityQueue {
   public:
    using key_type = Key;
    using mapped_type = T;
    using value_type = std::pair<key_type, mapped_type>;
    using key_compare = std::conditional_t<Min, std::greater<key_type>, std::less<key_type>>;

    class value_compare {
        friend TBBPriorityQueue;
        [[no_unique_address]] key_compare comp;

        explicit value_compare(key_compare const& c) : comp{c} {
        }

       public:
        constexpr bool operator()(value_type const& lhs, value_type const& rhs) const noexcept {
            return comp(lhs.first, rhs.first);
        }
    };

   private:
    using pq_type = tbb::concurrent_priority_queue<value_type, value_compare>;

   public:
    using handle_type = TBBPriorityQueue&;

   private:
    pq_type pq_;

   public:
    TBBPriorityQueue(std::size_t initial_capacity) : pq_(initial_capacity, value_compare{key_compare{}}) {
    }

    void push(value_type const& value) {
        pq_->push(value);
    }

    std::optional<value_type> try_pop() {
        value_type retval;
        if (!pq_->try_pop(retval)) {
            return std::nullopt;
        }
        return retval;
    }

    handle_type get_handle() {
        return *this;
    }
};

template <bool Min = true, typename Key = unsigned long, typename T = unsigned long>
using PQWrapper = TBBPriorityQueue<Min, Key, T>;

inline void add_options(cxxopts::Options& /*options*/) {
}

template <bool Min = true, typename Key = unsigned long, typename T = unsigned long>
TBBPriorityQueue<Min, Key, T> create(int /*num_threads*/, std::size_t initial_capacity, cxxopts::ParseResult const& /*result*/) {
    return TBBPriorityQueue<Min, Key, T>{initial_capacity};
}

template <typename PQ>
std::ostream& describe(PQ const& /*unused*/, std::ostream& out) {
    out << "TBB PriorityQueue";
    return out;
}
}  // namespace wrapper
