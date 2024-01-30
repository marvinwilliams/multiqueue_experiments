#pragma once

#include "util.hpp"

#include <cxxopts.hpp>

#include <tbb/concurrent_priority_queue.h>

#include <optional>
#include <ostream>
#include <utility>

namespace wrapper::tbb_pq {

template <bool Min, typename Key = unsigned long, typename T = unsigned long>
class TBBPriorityQueue {
   public:
    using key_type = Key;
    using mapped_type = T;
    using value_type = std::pair<key_type, mapped_type>;
    using key_compare = std::conditional_t<Min, std::greater<key_type>, std::less<key_type>>;
    using value_compare = util::ValueCompare<value_type, key_compare, util::PairFirst>;

   private:
    using pq_type = tbb::concurrent_priority_queue<value_type, value_compare>;

    pq_type pq_;

   public:
    using handle_type = util::SelfHandle<TBBPriorityQueue>;
    using settings_type = util::EmptySettings;
    TBBPriorityQueue(int /*unused*/, std::size_t initial_capacity, settings_type const& /*unused*/)
        : pq_(initial_capacity, value_compare{key_compare{}}) {
    }

    void push(value_type const& value) {
        pq_.push(value);
    }

    std::optional<value_type> try_pop() {
        value_type retval;
        if (!pq_.try_pop(retval)) {
            return std::nullopt;
        }
        return retval;
    }

    static void write_human_readable(std::ostream& out) {
        out << "TBB Priority Queue";
    }

    handle_type get_handle() {
        return handle_type{*this};
    }
};

}  // namespace wrapper::tbb_pq
