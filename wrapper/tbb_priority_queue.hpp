#pragma once
#ifndef WRAPPER_TBB_PRIORITY_QUEUE_HPP_INCLUDED
#define WRAPPER_TBB_PRIORITY_QUEUE_HPP_INCLUDED

#include <tbb/concurrent_priority_queue.h>

#include <string>
#include <utility>

namespace wrapper {

template <typename KeyType, typename T>
class TBBPriorityQueue {
   public:
    using key_type = KeyType;
    using mapped_type = T;
    using value_type = std::pair<key_type, mapped_type>;
    struct value_compare {
        bool operator()(value_type const& lhs,
                        value_type const& rhs) const noexcept {
            return lhs.first > rhs.first;
        }
    };

   private:
    using pq_type = tbb::concurrent_priority_queue<value_type, value_compare>;

   public:
    class Handle {
        friend TBBPriorityQueue;
        pq_type* pq_;

       public:
        void push(value_type const& value) { pq_->push(value); }
        bool try_pop(value_type& retval) {
            return pq_->try_pop(retval);
        }
    };

   private:
    pq_type pq_;

   public:
    TBBPriorityQueue(std::size_t capacity, unsigned int /* num_threads */)
        : pq_(capacity) {}

    Handle get_handle() {
        auto h = Handle{};
        h.pq_ = &pq_;
        return h;
    }

    void push(value_type const& value) { pq_.push(value); }
    bool try_pop(value_type& retval) { return pq_.try_pop(retval); }

    static std::string description() { return "TBBPriorityQueue"; }
};

}  // namespace wrapper

#endif
