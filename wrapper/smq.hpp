#pragma once

#include "StealingMultiQueue.hpp"
#include "util.hpp"


#include <atomic>
#include <cstdio>
#include <optional>
#include <ostream>
#include <utility>

namespace wrapper::stealing_mq {

template <bool Min, typename Key = unsigned long, typename T = Key>
class StealingMQ {
   public:
    using key_type = Key;
    using mapped_type = T;
    using value_type = std::pair<key_type, mapped_type>;
    using key_compare = std::conditional_t<Min, std::greater<>, std::less<>>;
    using value_compare = util::ValueCompare<value_type, key_compare, util::PairFirst>;

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
        bool push(value_type const& value) {
            value_type v_arr[1] = {value};
            pq_->push(id_, v_arr, v_arr + 1);
            return true;
        }

        std::optional<value_type> try_pop() {
            return pq_->pop(id_);
        }
    };

   public:
    using handle_type = Handle;
    using settings_type = util::EmptySettings;

   private:
    pq_type pq_;

   public:
    explicit StealingMQ(int num_threads, std::size_t /*unused*/, settings_type const& /*unused*/) : pq_(num_threads) {
    }

    static void write_human_readable(std::ostream& out) {
        out << "StealingMQ\n";
        out << "  Steal probability: " << StealProb << '\n';
        out << "  Steal batch size: " << StealBatchSize << '\n';
    }

    handle_type get_handle() {
        static std::atomic_int id{0};
        return Handle{&pq_, id.fetch_add(1)};
    }
};

}  // namespace wrapper::stealing_mq
