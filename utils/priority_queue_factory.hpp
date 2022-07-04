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
#ifndef UTILS_PRIORITY_QUEUE_FACTORY_HPP_INCLUDED
#define UTILS_PRIORITY_QUEUE_FACTORY_HPP_INCLUDED

#if defined PQ_MQ
#include "multiqueue/multiqueue.hpp"
#elif defined PQ_MF
#include "multififo/multififo.hpp"
#elif defined PQ_CAPQ
#include "wrapper/capq.hpp"
#elif defined PQ_KLSM256 || defined PQ_KLSM1024 || defined PQ_KLSM4096
#include "wrapper/klsm.hpp"
#elif defined PQ_LINDEN
#include "wrapper/linden.hpp"
#elif defined PQ_SPRAYLIST
#include "wrapper/spraylist.hpp"
#elif defined PQ_TBB_Q
#include "wrapper/tbb_queue.hpp"
#elif defined PQ_TBB_PQ
#include "wrapper/tbb_priority_queue.hpp"
#else
#error No valid PQ specified
#endif

#include <cstdint>
#include <functional>
#include <optional>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>

namespace util {
namespace detail {

template <typename KeyType, typename ValueType>
struct PriorityQueueTypeFactory {
#if defined PQ_MQ
    using type = multiqueue::MultiQueue<KeyType, ValueType>;
#elif defined PQ_MF
    using type = multififo::MultiFifo<std::pair<KeyType, ValueType>>;
#elif defined PQ_CAPQ
#error Not available with generic types
#elif defined PQ_KLSM256
    using type = wrapper::Klsm<KeyType, ValueType, 256>;
#elif defined PQ_KLSM1024
    using type = wrapper::Klsm<KeyType, ValueType, 1024>;
#elif defined PQ_KLSM4096
    using type = wrapper::Klsm<KeyType, ValueType, 4096>;
#elif defined PQ_LINDEN
#error Not available with generic types
#elif defined PQ_SPRAYLIST
#error Not available with generic types
#elif defined PQ_TBB_Q
    using type = wrapper::TBBQueue<KeyType, ValueType>;
#elif defined PQ_TBB_PQ
    using type = wrapper::TBBPriorityQueue<KeyType, ValueType>;
#endif
};

template <>
struct PriorityQueueTypeFactory<unsigned long, unsigned long> {
    using KeyType = unsigned long;
    using ValueType = unsigned long;
#if defined PQ_MQ
    using type = multiqueue::MultiQueue<KeyType, ValueType>;
#elif defined PQ_MF
    using type = multififo::MultiFifo<std::pair<unsigned long, unsigned long>>;
#elif defined PQ_CAPQ
    using type = wrapper::Capq<true, true, true>;
#elif defined PQ_KLSM256
    using type = wrapper::Klsm<KeyType, ValueType, 256>;
#elif defined PQ_KLSM1024
    using type = wrapper::Klsm<KeyType, ValueType, 1024>;
#elif defined PQ_KLSM4096
    using type = wrapper::Klsm<KeyType, ValueType, 4096>;
#elif defined PQ_LINDEN
    using type = wrapper::Linden;
#elif defined PQ_SPRAYLIST
    using type = wrapper::Spraylist;
#elif defined PQ_TBB_Q
    using type = wrapper::TBBQueue<KeyType, ValueType>;
#elif defined PQ_TBB_PQ
    using type = wrapper::TBBPriorityQueue<KeyType, ValueType>;
#endif
};

template <typename PriorityQueue, typename = void>
struct has_params : std::false_type {};

template <typename PriorityQueue>
struct has_params<PriorityQueue,
                  std::void_t<typename PriorityQueue::config_type>>
    : std::true_type {};

template <typename PriorityQueue>
static constexpr bool has_params_v = has_params<PriorityQueue>::value;

}  // namespace detail

struct PriorityQueueParameters {
    std::optional<std::uint64_t> seed;
    std::optional<std::size_t> c;
    std::optional<unsigned int> stickiness;
};

template <typename KeyType, typename ValueType>
struct PriorityQueueFactory {
    using type =
        typename detail::PriorityQueueTypeFactory<KeyType, ValueType>::type;
};

template <typename PriorityQueue>
static PriorityQueue create_pq(std::size_t initial_capacity,
                               unsigned int num_threads,
                               PriorityQueueParameters const& params = {}) {
    if constexpr (detail::has_params_v<PriorityQueue>) {
        typename PriorityQueue::config_type config{};
        if (params.seed) {
            config.seed = *params.seed;
        }
        if (params.c) {
            config.c = *params.c;
        }
#ifdef MQ_HAS_STICKINESS
        if (params.stickiness) {
            config.stick_policy_config.stickiness = *params.stickiness;
        }
#endif
        return PriorityQueue(initial_capacity, num_threads, config);
    } else {
        return PriorityQueue(initial_capacity, num_threads);
    }
}

}  // namespace util

#endif  //! UTILS_PRIORITY_QUEUE_FACTORY_HPP_INCLUDED
