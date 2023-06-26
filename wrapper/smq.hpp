#pragma once

#include "StealingMultiQueue.hpp"

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

template <typename Key, typename Value, typename Compare = std::less<Key>, std::size_t StealProb = 8, std::size_t StealBatchSize = 8>
class StealingMQ {
   public:
    using key_type = Key;
    using value_type = Value;
    struct config_type {};
    using key_comp = Compare;

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
        pq_type& pq_;
        int id_;

       public:
        Handle(pq_type& pq, int id) : pq_{pq}, id_{id} {
        }

        void push(value_type const& value) {
            pq_.push(id_, &value, (&value) + 1);
        }

        std::optional<value_type> try_pop() {
            return pq_.pop(id_);
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
        return Handle{*pq_.get(), id.fetch_add(1)};
    }

    std::ostream& describe(std::ostream& out) {
        out << "Stealing MultiQueue (StealProb=" << StealProb << ", StealBatchSize=" << StealBatchSize << ")";
        return out;
    }
};

}  // namespace wrapper
