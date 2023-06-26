#pragma once

#include "cxxopts.hpp"

#include <tbb/concurrent_queue.h>

#include <optional>
#include <ostream>
#include <utility>

namespace wrapper {

// Min is ignored since this is a FIFO
template <typename KeyType, typename Value>
class TBBQueue {
   public:
    using key_type = KeyType;
    using value_type = Value;
    struct config_type {};

   private:
    using pq_type = tbb::concurrent_queue<value_type>;

   public:
    class Handle {
        friend TBBQueue;
        pq_type* pq_;

       public:
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
    };

    using handle_type = Handle;

   private:
    pq_type pq_;

   public:
    static void add_options(cxxopts::Options& /*options*/, config_type& /*config*/) {
    }

    TBBQueue(int /*num_threads*/, std::size_t /*initial_capacity*/, config_type const& /*options*/) {
    }

    Handle get_handle() {
        auto h = Handle{};
        h.pq_ = &pq_;
        return h;
    }

    void push(value_type const& value) {
        pq_.push(value);
    }

    bool try_pop(value_type& retval) {
        return pq_.try_pop(retval);
    }

    std::ostream& describe(std::ostream& out) {
        out << "TBB Queue";
        return out;
    }
};

}  // namespace wrapper
