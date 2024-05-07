#pragma once

#include "util.hpp"
#define MQ_MODE_RANDOM
#include "wrapper/multiqueue.hpp"
#undef MQ_MODE_RANDOM

#include <mutex>

#ifdef __SSE2__
#include <emmintrin.h>
#define PAUSE _mm_pause()
#else
#define PAUSE void(0)
#endif
namespace wrapper::mq_pq {

template <bool Min, typename Key = unsigned long, typename T = Key>
class MultiQueuePQ {
   public:
    using key_type = Key;
    using mapped_type = T;
    using value_type = std::pair<key_type, mapped_type>;
    using key_compare = std::conditional_t<Min, std::greater<>, std::less<>>;
    using value_compare = util::ValueCompare<value_type, key_compare, util::PairFirst>;

   private:
    using pq_type = typename wrapper::multiqueue::MultiQueue<Min, Key, T>::pq_type;

    std::mutex lock_{};
    pq_type pq_{};

   public:
    using handle_type = util::SelfHandle<MultiQueuePQ>;
    using settings_type = util::EmptySettings;

    MultiQueuePQ(int /*unused*/, std::size_t initial_capacity, settings_type const& /*unused*/) {
        pq_.reserve(initial_capacity);
    }

    void push(value_type const& value) {
        std::scoped_lock<std::mutex> guard(lock_);
        pq_.push(value);
    }

    std::optional<value_type> try_pop() {
        std::scoped_lock<std::mutex> guard(lock_);
        if (pq_.empty()) {
            return std::nullopt;
        }
        auto retval = pq_.top();
        pq_.pop();
        return retval;
    }

    static void write_human_readable(std::ostream& out) {
        out << "MultiQueuePQ" << '\n';
#ifdef MQ_USE_BTREE
        out << "  PQ: tlx::btree" << '\n';
#else
#ifdef MQ_USE_STD_PQ
        out << "  PQ: std::priority_queue" << '\n';
#else
        out << "  PQ: d-ary heap" << '\n';
        out << "  Heap arity: " << wrapper::multiqueue::heap_arity << '\n';
#endif
        out << "  Insertion buffer size: " << wrapper::multiqueue::insertion_buffer_size << '\n';
        out << "  Deletion buffer size: " << wrapper::multiqueue::deletion_buffer_size << '\n';
#endif
    }

    handle_type get_handle() {
        return handle_type{*this};
    }
};

}  // namespace wrapper::mq_pq

#undef PAUSE
