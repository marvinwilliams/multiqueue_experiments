#pragma once

#include "util.hpp"
#define MQ_MODE_RANDOM
#include "wrapper/multiqueue.hpp"
#undef MQ_MODE_RANDOM

#include <atomic>

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

    std::atomic_bool lock_{};
    pq_type pq_{};

   public:
    using handle_type = util::SelfHandle<MultiQueuePQ>;
    using settings_type = util::EmptySettings;

   private:
    void lock() {
        /* while (true) { */
        /*     while (lock_.load(std::memory_order_relaxed)) { */
        /*         PAUSE; */
        /*     } */
        /*     if (!lock_.exchange(true, std::memory_order_acquire)) { */
        /*         return; */
        /*     } */
        /* } */
        while (true) {
            if (!lock_.exchange(true, std::memory_order_acquire)) {
                return;
            }
            do {
                PAUSE;
            } while (lock_.load(std::memory_order_relaxed));
        }
    }

    void unlock() {
        lock_.store(false, std::memory_order_release);
    }

   public:
    MultiQueuePQ(int /*unused*/, std::size_t initial_capacity, settings_type const& /*unused*/) {
        pq_.reserve(initial_capacity);
    }

    void push(value_type const& value) {
        lock();
        pq_.push(value);
        unlock();
    }

    std::optional<value_type> try_pop() {
        lock();
        if (pq_.empty()) {
            unlock();
            return std::nullopt;
        }
        auto retval = pq_.top();
        pq_.pop();
        unlock();
        return retval;
    }

    static void write_human_readable(std::ostream& out) {
        out << "MultiQueue PQ\n";
    }

    handle_type get_handle() {
        return handle_type{*this};
    }
};

}  // namespace wrapper::mq_pq

#undef PAUSE
