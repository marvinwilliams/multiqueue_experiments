#pragma once

#include "StealingMultiQueue.hpp"

#include "wrapper/priority.hpp"

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

namespace wrapper {

template <typename Key, typename T, Priority P = Priority::Min>
class StealingMQ {
   public:
    using key_type = Key;
    using mapped_type = T;
    using value_type = std::pair<key_type, mapped_type>;
    struct config_type {};
    using key_comp = std::conditional_t<P == Priority::Min, std::greater<Key>, std::less<Key>>;

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
    class value_compare {
        [[no_unique_address]] key_comp comp;

       public:
        explicit value_compare(key_comp const& compare = key_comp{}) : comp{compare} {
        }

        constexpr bool operator()(value_type const& lhs, value_type const& rhs) const noexcept {
            return comp(lhs.first, rhs.first);
        }
    };
    using pq_type = Galois::WorkList::StealingMultiQueue<value_type, value_compare, StealProb, StealBatchSize>;

   public:
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

    using handle_type = Handle;

   private:
    alignas(64) std::unique_ptr<pq_type> pq_;

   public:
    static void add_options(cxxopts::Options& /*options*/, config_type& /*config*/) {
    }

    explicit StealingMQ(int num_threads, std::size_t /*initial_capacity*/, config_type const& /*config*/)
        : pq_(new pq_type(num_threads)) {
    }

    Handle get_handle() {
        static std::atomic_int id{0};
        return Handle{pq_.get(), id.fetch_add(1)};
    }

    std::ostream& describe(std::ostream& out) {
        out << "Stealing MultiQueue (StealProb=" << StealProb << ", StealBatchSize=" << StealBatchSize << ")";
        return out;
    }
};

}  // namespace wrapper
