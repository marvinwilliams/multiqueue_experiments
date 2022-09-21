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
#include "multiqueue/buffered_pq.hpp"
#include "multiqueue/build_config.hpp"
#include "multiqueue/config.hpp"
#include "multiqueue/heap.hpp"
#include "multiqueue/multiqueue.hpp"
#include "multiqueue/stick_policy.hpp"
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

#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <ostream>
#include <queue>
#include <type_traits>
#include <utility>

namespace util {

struct PriorityQueueConfig {
    std::optional<std::uint64_t> seed;
    std::optional<int> c;
    std::optional<int> stickiness;
};

namespace detail {

#if defined PQ_MQ
static constexpr ::multiqueue::StickPolicy stick_policy =
#ifdef MQ_STICK_POLICY_NONE
    multiqueue::StickPolicy::None
#elif defined MQ_STICK_POLICY_RANDOM
    multiqueue::StickPolicy::Random
#elif defined MQ_STICK_POLICY_RANDOM_STRICT
    multiqueue::StickPolicy::RandomStrict
#elif defined MQ_STICK_POLICY_SWAPPING
    multiqueue::StickPolicy::Swapping
#elif defined MQ_STICK_POLICY_SWAPPING_LAZY
    multiqueue::StickPolicy::SwappingLazy
#elif defined MQ_STICK_POLICY_SWAPPING_BLOCKING
    multiqueue::StickPolicy::SwappingBlocking
#elif defined MQ_STICK_POLICY_PERMUTATION
    multiqueue::StickPolicy::Permutation
#else
    multiqueue::StickPolicy::None
#endif
    ;

template <typename T, typename Compare>
using InnerPriorityQueue =
#ifndef MQ_DISABLE_BUFFERING
    multiqueue::BufferedPQ<
#endif
#ifdef MQ_PQ_STD
        std::priority_queue<T, std::vector<T>, Compare>
#else
    multiqueue::Heap<T, Compare
#ifdef MQ_HEAP_ARITY
                     ,
                     MQ_HEAP_ARITY
#endif
                     >
#endif
#ifndef MQ_DISABLE_BUFFERING
        ,
#ifdef MQ_INSERTION_BUFFERSIZE
        MQ_INSERTION_BUFFERSIZE
#else
        multiqueue::BuildConfiguration::DefaultInsertionBuffersize
#endif
        ,
#ifdef MQ_DELETION_BUFFERSIZE
        MQ_DELETION_BUFFERSIZE
#else
        multiqueue::BuildConfiguration::DefaultDeletionBuffersize
#endif
        >
#endif
    ;
#endif

template <typename KeyType, typename ValueType, bool Min = true>
struct PriorityQueueTypeFactory {
#if defined PQ_MQ
    using type = multiqueue::MultiQueue<KeyType, ValueType, std::greater<KeyType>, stick_policy, InnerPriorityQueue>;
#elif defined PQ_MF
    using type = multififo::MultiFifo<std::pair<KeyType, ValueType>>;
#elif defined PQ_KLSM256
    using type = wrapper::Klsm<KeyType, ValueType, 256>;
#elif defined PQ_KLSM1024
    using type = wrapper::Klsm<KeyType, ValueType, 1024>;
#elif defined PQ_KLSM4096
    using type = wrapper::Klsm<KeyType, ValueType, 4096>;
#elif defined PQ_TBB_Q
    using type = wrapper::TBBQueue<KeyType, ValueType>;
#elif defined PQ_TBB_PQ
    using type = wrapper::TBBPriorityQueue<KeyType, ValueType, std::greater<>>;
#endif
};

template <typename KeyType, typename ValueType>
struct PriorityQueueTypeFactory<KeyType, ValueType, false> {
#if defined PQ_MQ
    using type = multiqueue::MultiQueue<KeyType, ValueType, std::less<KeyType>, stick_policy, InnerPriorityQueue>;
#elif defined PQ_MF
    using type = multififo::MultiFifo<std::pair<KeyType, ValueType>>;
#elif defined PQ_TBB_Q
    using type = wrapper::TBBQueue<KeyType, ValueType>;
#elif defined PQ_TBB_PQ
    using type = wrapper::TBBPriorityQueue<KeyType, ValueType>;
#endif
};

template <>
struct PriorityQueueTypeFactory<unsigned long, unsigned long, true> {
    using KeyType = unsigned long;
    using ValueType = unsigned long;
#if defined PQ_MQ
    using type = multiqueue::MultiQueue<KeyType, ValueType, std::greater<>, stick_policy, InnerPriorityQueue>;
#elif defined PQ_MF
    using type = multififo::MultiFifo<std::pair<KeyType, ValueType>>;
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
    using type = wrapper::TBBPriorityQueue<KeyType, ValueType, std::greater<>>;
#endif
};

template <>
struct PriorityQueueTypeFactory<unsigned long, unsigned long, false> {
    using KeyType = unsigned long;
    using ValueType = unsigned long;
#if defined PQ_MQ
    using type = multiqueue::MultiQueue<KeyType, ValueType, std::less<>, stick_policy, InnerPriorityQueue>;
#elif defined PQ_MF
    using type = multififo::MultiFifo<std::pair<KeyType, ValueType>>;
#elif defined PQ_TBB_Q
    using type = wrapper::TBBQueue<KeyType, ValueType>;
#elif defined PQ_TBB_PQ
    using type = wrapper::TBBPriorityQueue<KeyType, ValueType>;
#endif
};

template <typename T>
struct DescribeType {};

template <typename T, typename Container, typename Compare>
std::ostream &describe_type(std::ostream &out, DescribeType<std::priority_queue<T, Container, Compare>> /*unused*/) {
    out << "std::priority_queue\n";
    return out;
}

#ifdef PQ_MQ
template <typename PQ, std::size_t I, std::size_t D>
std::ostream &describe_type(std::ostream &out, DescribeType<multiqueue::BufferedPQ<PQ, I, D>> /*unused*/) {
    out << "Buffered PQ\n";
    out << "Buffer sizes (Insertion/Deletion): " << I << '/' << D << '\n';
    out << "Inner PQ type: ";
    describe_type(out, DescribeType<PQ>{});
    return out;
}

template <typename T, typename Compare, unsigned int Degree, typename Container>
std::ostream &describe_type(std::ostream &out,
                            DescribeType<multiqueue::Heap<T, Compare, Degree, Container>> /*unused*/) {
    out << "Heap\n";
    out << "Degree: " << Degree << '\n';
    return out;
}
#endif
}  // namespace detail

template <typename KeyType, typename ValueType, bool Min>
struct PriorityQueueFactory {
    using type = typename detail::PriorityQueueTypeFactory<KeyType, ValueType, Min>::type;
    static type create(unsigned int num_threads, PriorityQueueConfig const &params = {}) {
#ifdef PQ_MQ
        multiqueue::Config config{};
        if (params.seed) {
            config.seed = *params.seed;
        }
        if (params.c) {
            config.c = *params.c;
        }
        if (params.stickiness) {
            config.stickiness = *params.stickiness;
        }
        return type(num_threads, config);
#else
        return type(num_threads);
#endif
    }

    static std::ostream &describe(std::ostream &out, type const &pq) {
        pq.describe(out);
#ifdef PQ_MQ
        out << "Inner pq: ";
        detail::describe_type(out, detail::DescribeType<typename type::inner_pq_type>{});
#endif
        return out;
    }
};

}  // namespace util
