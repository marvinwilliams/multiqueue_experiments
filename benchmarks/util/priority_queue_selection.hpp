/**
******************************************************************************
* @file:   select_queue.hpp
*
* @author: Marvin Williams
* @date:   2021/02/22 13:23
* @brief:
*******************************************************************************
**/
#pragma once

#if defined PQ_MQ
#include "wrapper/multiqueue.hpp"
template <typename Key, typename T, bool Min>
using PriorityQueue = wrapper::MultiQueue<Key, T, Min>;
#elif defined PQ_CAPQ
#include "wrapper/capq.hpp"
template <typename Key, typename T, bool Min>
using PriorityQueue = wrapper::CAPQ<Key, T, Min>;
#elif defined PQ_KLSM
#include "wrapper/klsm.hpp"
template <typename Key, typename T, bool Min>
using PriorityQueue = wrapper::KLsm<Key, T, Min>;
#elif defined PQ_LINDEN
#include "wrapper/linden.hpp"
template <typename Key, typename T, bool Min>
using PriorityQueue = wrapper::Linden<Key, T, Min>;
#elif defined PQ_SPRAYLIST
#include "wrapper/spraylist.hpp"
template <typename Key, typename T, bool Min>
using PriorityQueue = wrapper::Spraylist<Key, T, Min>;
#elif defined PQ_TBB_PQ
#include "wrapper/tbb_priority_queue.hpp"
template <typename Key, typename T, bool Min>
using PriorityQueue = wrapper::TBBPriorityQueue<Key, T, Min>;
#elif defined PQ_TBB_FIFO
#include "wrapper/tbb_queue.hpp"
template <typename Key, typename T, bool Min>
using PriorityQueue = wrapper::TBBQueue<Key, T, Min>;
#elif defined PQ_SMQ
#include "wrapper/smq.hpp"
template <typename Key, typename T, bool Min>
using PriorityQueue = wrapper::StealingMultiQueue<Key, T, Min>;
#else
#error No valid PQ specified
#endif

#include "operation_log.hpp"

#include <chrono>
#include <ctime>

struct OperationCounters {
    long long pop = 0;
    long long push = 0;
    long long failed_pop = 0;
};

template <typename PriorityQueue, bool LogOperations /*= false*/>
class Handle {
   public:
    using key_type = typename PriorityQueue::key_type;
    using value_type = typename PriorityQueue::value_type;

   private:
    typename PriorityQueue::handle_type handle_;
    OperationCounters counters_;

   public:
    explicit Handle(int /*id*/, PriorityQueue& pq) : handle_(pq.get_handle()) {
    }

    void push(key_type const& key) {
        handle_.push({key, key});
        ++counters_.push;
    }

    std::optional<value_type> try_pop() {
        auto ret = handle_.try_pop();
        if (!ret) {
            ++counters_.failed_pop;
            return std::nullopt;
        }
        ++counters_.pop;
        return *ret;
    }

    auto const& get_counters() const {
        return counters_;
    }

    auto const& get_stats() const {
        return this->handle_.get_stats();
    }
};

template <typename PriorityQueue>
class Handle<PriorityQueue, true> {
   public:
    using key_type = typename PriorityQueue::key_type;
    using value_type = typename PriorityQueue::value_type;

   private:
    typename PriorityQueue::handle_type handle_;
    unsigned long value_;
    operation_log::OperationLog log_;

    std::uint64_t get_tick() noexcept {
        timespec ts{};
        clock_gettime(CLOCK_MONOTONIC_RAW, &ts);
        return static_cast<std::uint64_t>(
            (std::chrono::seconds(ts.tv_sec) + std::chrono::nanoseconds(ts.tv_nsec)).count());
    }

   public:
    explicit Handle(int id, PriorityQueue& pq) : handle_(pq.get_handle()), value_(operation_log::pack(id, 0)) {
    }

    void reserve_push_log(std::size_t size) {
        log_.pushes.reserve(size);
    }

    void reserve_pop_log(std::size_t size) {
        log_.pops.reserve(size);
    }

    void push(key_type const& key) {
        handle_.push({key, value_});
        auto tick = get_tick();
        log_.pushes.push_back({tick, key});
        ++value_;
    }

    std::optional<value_type> try_pop() {
        auto tick = get_tick();
        auto retval = handle_.try_pop();
        if (!retval) {
            return std::nullopt;
        }
        log_.pops.push_back({tick, retval->second});
        return *retval;
    }

    operation_log::OperationLog const& log() const& noexcept {
        return log_;
    }

    auto const& get_stats() const {
        return this->handle_.get_stats();
    }
};
