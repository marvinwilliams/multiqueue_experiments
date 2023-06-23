#pragma once

#include <algorithm>
#include <chrono>

#include "Galois/WorkList/StealingMultiQueue.h"

#include <algorithm>
#include <cstdio>
#include <limits>
#include <memory>
#include <optional>
#include <ostream>
#include <utility>

namespace wrapper {

template <typename Key, typename T, typename Compare = std::less<Key>, std::size_t StealProb = 8,
          std::size_t StealBatchSize = 8>
class StealingMQ {
   public:
    using key_type = Key;
    using mapped_type = T;
    using value_type = std::pair<key_type, mapped_type>;

   private:
    class value_compare {
        [[no_unique_address]] Compare comp;

       public:
        explicit value_compare(Compare const& compare = Compare{}) : comp{compare} {
        }

        constexpr bool operator()(value_type const& lhs, value_type const& rhs) const noexcept {
            return comp(lhs.first, rhs.first);
        }
    };
    using pq_type = Galois::WorkList::StealingMultiQueue<value_type, value_compare, StealProb, StealBatchSize>;

   public:
    class Handle {
        pq_type& pq_;
        int id_;

       public:
        Handle(pq_type& pq, int id) : pq_{pq}, id_{id} {
        }

        void push(value_type const& value) {
            pq_.push(id_, &value, (&value) + 1);
        }

        template <typename Iter>
        unsigned int push(Iter b, Iter e) {
            return pq_.push(id_, b, e);
        }

        std::optional<value_type> try_pop() {
            auto t = pq_.pop(id_);
            if (!t.is_initialized()) {
                return std::nullopt;
            }
            return t.get();
        }
    };

    using handle_type = Handle;

   private:
    alignas(64) std::unique_ptr<pq_type> pq_;

   public:
    StealingMQ(int num_threads, std::size_t initial_capacity);

    Handle get_handle(int id);
};

template <typename Key, typename T, typename Compare, std::size_t StealProb, std::size_t StealBatchSize>
StealingMQ<Key, T, Compare, StealProb, StealBatchSize>::StealingMQ(int num_threads, std::size_t /*unused*/)
    : pq_(new pq_type(num_threads)) {
}

template <typename Key, typename T, typename Compare, std::size_t StealProb, std::size_t StealBatchSize>
typename StealingMQ<Key, T, Compare, StealProb, StealBatchSize>::Handle
StealingMQ<Key, T, Compare, StealProb, StealBatchSize>::get_handle(int id) {
    return Handle{*pq_.get(), id};
}
}  // namespace wrapper
