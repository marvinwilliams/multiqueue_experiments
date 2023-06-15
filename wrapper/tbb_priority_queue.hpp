#pragma once

#include "cxxopts.hpp"

#include <tbb/concurrent_priority_queue.h>

#include <ostream>
#include <utility>

namespace wrapper {

template <typename KeyType, typename T, bool Min>
class TBBPriorityQueue {
   public:
    using key_type = KeyType;
    using mapped_type = T;
    using value_type = std::pair<key_type, mapped_type>;
    using key_compare = std::conditional_t<Min, std::greater<key_type>, std::less<key_type>>;
    class value_compare {
        friend class TBBPriorityQueue<KeyType, T, key_compare>;
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
    class Handle {
        friend TBBPriorityQueue;
        pq_type* pq_;

       public:
        void push(value_type const& value) {
            pq_->push(value);
        }
        bool try_pop(value_type& retval) {
            return pq_->try_pop(retval);
        }
    };

    using handle_type = Handle;

   private:
    pq_type pq_;

   public:
    static void add_options(cxxopts::Options& /*options*/) {
    }

    TBBPriorityQueue(int /*num_threads*/, std::size_t initial_capacity, cxxopts::ParseResult const& /*options*/)
        : pq_(initial_capacity, value_compare{key_compare{}}) {
    }

    Handle get_handle() {
        auto h = Handle{};
        h.pq_ = &pq_;
        return h;
    }

    std::ostream& describe(std::ostream& out) {
        out << "TBB PriorityQueue";
        return out;
    }
};

}  // namespace wrapper
