#pragma once

#include "StealingMultiQueue.hpp"
#include "util.hpp"

#include "cxxopts.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <limits>
#include <memory>
#include <optional>
#include <ostream>
#include <utility>

namespace wrapper::stealing_mq {

template <bool Min = true, typename Key = unsigned long, typename Value = std::pair<unsigned long, unsigned long>,
          typename KeyOfValue = util::KeyOfValue<Key, Value>>
class StealingMQ {
   public:
    using key_type = Key;
    using value_type = Value;
    using key_compare = std::conditional_t<Min, std::greater<key_type>, std::less<key_type>>;
    using value_compare = util::ValueCompare<value_type, key_compare, KeyOfValue>;

#ifdef SMQ_STEAL_PROB
    static constexpr std::size_t StealProb = SMQ_STEAL_PROB;
#else
    static constexpr std::size_t StealProb = 8;
#endif
#ifdef SMQ_STEAL_BATCH_SIZE
    static constexpr std::size_t StealBatchSize = SMQ_STEAL_BATCH_SIZE;
#else
    static constexpr std::size_t StealBatchSize = 64;
#endif

   private:
    using pq_type = Galois::WorkList::StealingMultiQueue<value_type, value_compare, StealProb, StealBatchSize>;

    class Handle {
        friend StealingMQ;

        pq_type* pq_;
        int id_;

        Handle(pq_type* pq, int id) : pq_{pq}, id_{id} {
        }

       public:
        void push(value_type const& value) {
            pq_->push(id_, &value, (&value) + 1);
        }

        std::optional<value_type> try_pop() {
            return pq_->pop(id_);
        }
    };

   public:
    using handle_type = Handle;

   private:
    pq_type pq_;

   public:
    explicit StealingMQ(int num_threads) : pq_(num_threads) {
    }

    Handle get_handle() {
        static std::atomic_int id{0};
        return Handle{&pq_, id.fetch_add(1)};
    }
};

template <bool Min = true, typename Key = unsigned long, typename Value = std::pair<unsigned long, unsigned long>,
          typename KeyOfValue = util::KeyOfValue<Key, Value>>
using PQWrapper = StealingMQ<Min, Key, Value, KeyOfValue>;

inline void add_options(cxxopts::Options& /*options*/) {
}

template <bool Min = true, typename Key = unsigned long, typename Value = std::pair<unsigned long, unsigned long>,
          typename KeyOfValue = util::KeyOfValue<Key, Value>>
StealingMQ<Min, Key, Value, KeyOfValue> create(int num_threads, std::size_t /*initial_capacity*/,
                                               cxxopts::ParseResult const& /*result*/) {
    return StealingMQ<Min, Key, Value, KeyOfValue>{num_threads};
}

template <typename PQ>
std::ostream& describe(PQ const& pq, std::ostream& out) {
    out << "Stealing MultiQueue (StealProb=" << PQ::StealProb << ", StealBatchSize=" << PQ::StealBatchSize << ")";
    return out;
}

}  // namespace wrapper::stealing_mq
