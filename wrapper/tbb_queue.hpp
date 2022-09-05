#pragma once
#ifndef WRAPPER_TBB_QUEUE_HPP_INCLUDED
#define WRAPPER_TBB_QUEUE_HPP_INCLUDED

#include <tbb/concurrent_queue.h>

#include <string>
#include <utility>

namespace wrapper {

template <typename KeyType, typename T>
class TBBQueue {
   public:
    using key_type = KeyType;
    using mapped_type = T;
    using value_type = std::pair<key_type, mapped_type>;

   private:
    using pq_type = tbb::concurrent_queue<value_type>;

   public:
    class Handle {
        friend TBBQueue;
        pq_type* pq_;

       public:
        void push(value_type const& value) { pq_->push(value); }
        bool try_pop(value_type& retval) { return pq_->try_pop(retval); }
    };

    using handle_type = Handle;

   private:
    pq_type pq_;

   public:
    TBBQueue(unsigned int /* num_threads */) {}

    Handle get_handle() {
        auto h = Handle{};
        h.pq_ = &pq_;
        return h;
    }

    void push(value_type const& value) { pq_.push(value); }
    bool try_pop(value_type& retval) { return pq_.try_pop(retval); }

    static std::string description() { return "TBBQueue"; }
};

}  // namespace wrapper

#endif
