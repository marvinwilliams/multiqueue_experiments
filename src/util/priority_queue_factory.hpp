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
#include "multiqueue/heap.hpp"
#include "multiqueue/multiqueue.hpp"
#include "multiqueue/stick_policy.hpp"
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

#if defined PQ_MQ && defined USE_STD_PQ
#include <queue>
#endif

#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <ostream>
#include <type_traits>
#include <utility>

namespace util {

namespace detail {

#if defined PQ_MQ
static constexpr multiqueue::StickPolicy stick_policy =
#ifdef STICK_POLICY_NONE
    multiqueue::StickPolicy::None
#elif defined STICK_POLICY_RANDOM
    multiqueue::StickPolicy::Random
#elif defined STICK_POLICY_RANDOM_STRICT
    multiqueue::StickPolicy::RandomStrict
#elif defined STICK_POLICY_SWAPPING
    multiqueue::StickPolicy::Swapping
#elif defined STICK_POLICY_SWAPPING_LAZY
    multiqueue::StickPolicy::SwappingLazy
#elif defined STICK_POLICY_SWAPPING_BLOCKING
    multiqueue::StickPolicy::SwappingBlocking
#elif defined STICK_POLICY_PERMUTATION
    multiqueue::StickPolicy::Permutation
#else
    multiqueue::StickPolicy::None
#endif
    ;

template <typename T, typename Compare>
using SeqPriorityQueue = ::multiqueue::BufferedPQ<
#ifdef USE_STD_PQ
    std::priority_queue<T, std::vector<T>, Compare>
#elif defined HEAP_ARITY
    ::multiqueue::Heap<T, Compare, HEAP_ARITY>
#else
    ::multiqueue::Heap<T, Compare>
#endif
#if defined INSERTION_BUFFERSIZE && defined DELETION_BUFFERSIZE
    ,
    INSERTION_BUFFERSIZE, DELETION_BUFFERSIZE
#elif defined INSERTION_BUFFERSIZE || defined DELETION_BUFFERSIZE
#error Must specify either both buffersizes or none
#endif
    >;
#endif

template <typename KeyType, typename ValueType, bool Min = true>
struct PriorityQueueTypeFactory {
#if defined PQ_MQ
    using type = multiqueue::MultiQueue<KeyType, ValueType, std::greater<KeyType>, stick_policy, SeqPriorityQueue>;
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
    using type = multiqueue::MultiQueue<KeyType, ValueType, std::less<KeyType>, stick_policy, SeqPriorityQueue>;
#elif defined PQ_TBB_Q
    using type = wrapper::TBBQueue<KeyType, ValueType>;
#elif defined PQ_TBB_PQ
    using type = wrapper::TBBPriorityQueue<KeyType, ValueType, std::less<>>;
#endif
};

template <>
struct PriorityQueueTypeFactory<unsigned long, unsigned long, true> {
    using KeyType = unsigned long;
    using ValueType = unsigned long;
#if defined PQ_MQ
    using type = multiqueue::MultiQueue<KeyType, ValueType, std::greater<>, stick_policy, SeqPriorityQueue>;
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
    using type = multiqueue::MultiQueue<KeyType, ValueType, std::less<>, stick_policy, SeqPriorityQueue>;
#elif defined PQ_TBB_Q
    using type = wrapper::TBBQueue<KeyType, ValueType>;
#elif defined PQ_TBB_PQ
    using type = wrapper::TBBPriorityQueue<KeyType, ValueType, std::less<>>;
#endif
};

}  // namespace detail

template <typename KeyType, typename ValueType, bool Min>
using priority_queue_type = typename detail::PriorityQueueTypeFactory<KeyType, ValueType, Min>::type;

namespace describe {

template <typename T>
struct DescribeTag {};

template <typename T>
void describe_type(std::ostream &out, DescribeTag<T>) {
    out << "unknown\n";
}

#ifdef PQ_MQ
#ifdef USE_STD_PQ
template <typename T, typename Container, typename Compare>
void describe_type(std::ostream &out, DescribeTag<std::priority_queue<T, Container, Compare>>) {
    out << "std::priority_queue\n";
}
#endif

template <typename T, typename Compare, unsigned int Arity, typename Container>
void describe_type(std::ostream &out, DescribeTag<multiqueue::Heap<T, Compare, Arity, Container>>) {
    out << "Heap with arity " << Arity << '\n';
}

template <typename PQ, std::size_t I, std::size_t D>
void describe_type(std::ostream &out, DescribeTag<multiqueue::BufferedPQ<PQ, I, D>>) {
    out << "Buffered PQ\n";
    out << "Buffer sizes (Insertion/Deletion): " << I << '/' << D << '\n';
    out << "Backing PQ: ";
    describe_type(out, DescribeTag<PQ>{});
}

inline std::string stick_policy_to_name(multiqueue::StickPolicy policy) {
    switch (policy) {
        case multiqueue::StickPolicy::None:
            return "none";
        case multiqueue::StickPolicy::RandomStrict:
            return "random (strict)";
        case multiqueue::StickPolicy::Random:
            return "Random";
        case multiqueue::StickPolicy::Swapping:
            return "Swapping";
        case multiqueue::StickPolicy::SwappingBlocking:
            return "Swapping (blocking)";
        case multiqueue::StickPolicy::SwappingLazy:
            return "Swapping (lazy)";
        case multiqueue::StickPolicy::Permutation:
            return "Permutation";
        default:
            return "unknown";
    }
}

template <typename Key, typename T, typename KeyCompare, multiqueue::StickPolicy P,
          template <typename, typename> typename PriorityQueue, typename ValueTraits, typename SentinelTraits,
          typename Allocator>
void describe(
    std::ostream &out,
    multiqueue::MultiQueue<Key, T, KeyCompare, P, PriorityQueue, ValueTraits, SentinelTraits, Allocator> const &mq) {
    out << "MultiQueue\n";
    out << "Number of PQs: " << mq.num_pqs() << '\n';
    out << "Sentinel: " << SentinelTraits::sentinel() << " (" << (SentinelTraits::is_implicit ? "implicit" : "explicit")
        << ")\n";
    out << "Stick policy: " << stick_policy_to_name(P) << '\n';
    out << "Sequential PQ: ";
    describe_type(out,
                  DescribeTag<typename multiqueue::MultiQueue<Key, T, KeyCompare, P, PriorityQueue, ValueTraits,
                                                              SentinelTraits, Allocator>::pq_type::pq_type>{});
}

#else

template <typename PriorityQueue>
void describe(std::ostream &out, PriorityQueue const &pq) {
    pq.describe(out);
}

#endif
}  // namespace describe

}  // namespace util
