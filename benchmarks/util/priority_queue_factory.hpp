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
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#if defined PQ_MQ
#include "multiqueue/buffered_pq.hpp"
#include "multiqueue/config.hpp"
#include "multiqueue/multiqueue.hpp"
#include "multiqueue/stick_policy.hpp"

#include "cxxopts.hpp"

#include <iostream>

#define PQ_HAS_GENERIC
#define PQ_HAS_MAX
#define PQ_HAS_MAX_GENERIC
namespace detail {
#ifdef STICK_POLICY_NONE
static constexpr auto STICK_POLICY = ::multiqueue::StickPolicy::None;
#define STICK_POLICY_NAME none
#elif defined STICK_POLICY_RANDOM
static constexpr auto STICK_POLICY = ::multiqueue::StickPolicy::Random;
#define STICK_POLICY_NAME random
#elif defined STICK_POLICY_SWAPPING
static constexpr auto STICK_POLICY = ::multiqueue::StickPolicy::Swapping;
#define STICK_POLICY_NAME swapping
#elif defined STICK_POLICY_PERMUTATION
static constexpr auto STICK_POLICY = ::multiqueue::StickPolicy::Permutation;
#define STICK_POLICY_NAME permutation
#else
#error No valid stick policy specified
#endif
#ifdef USE_STD_PQ
#include <queue>
#else
#include "multiqueue/heap.hpp"
#endif
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
}  // namespace detail
template <typename Key, typename T>
using GenericMinPriorityQueue =
    multiqueue::MultiQueue<Key, T, std::greater<>, detail::STICK_POLICY, detail::SeqPriorityQueue>;
template <typename Key, typename T>
using GenericMaxPriorityQueue =
    multiqueue::MultiQueue<Key, T, std::less<>, detail::STICK_POLICY, detail::SeqPriorityQueue>;
using DefaultMinPriorityQueue = GenericMinPriorityQueue<unsigned long, unsigned long>;
using DefaultMaxPriorityQueue = GenericMaxPriorityQueue<unsigned long, unsigned long>;
static constexpr auto pq_name =
    "MultiQueue (buffer sizes: "
#if defined INSERTION_BUFFERSIZE && defined DELETION_BUFFERSIZE
    "i=" TOSTRING(INSERTION_BUFFERSIZE) " d=" TOSTRING(DELETION_BUFFERSIZE)
#else
                          "default"
#endif
        ", stick policy: " TOSTRING(STICK_POLICY_NAME)
        ", heap: "
#ifdef USE_STD_PQ
        "std"
#else
    "default"
#ifdef HEAP_ARITY
                                   " (d=" TOSTRING(HEAP_ARITY) ")"
#endif
#endif
        ")";
#undef STICK_POLICY_NAME

inline void add_pq_options(cxxopts::Options& options) {
    options.add_options()("c,factor", "The factor for queues", cxxopts::value<int>(), "NUMBER")
#ifndef PQ_MQ_NONE
        ("k,stickiness", "The stickiness period", cxxopts::value<int>(), "NUMBER");
#endif
}

using PriorityQueueConfig = multiqueue::Config;

inline PriorityQueueConfig get_pq_options(cxxopts::ParseResult const& result) {
    multiqueue::Config config;
    if (result.count("factor") > 0) {
        config.c = result["factor"].as<int>();
    }
#ifndef PQ_MQ_NONE
    if (result.count("stickiness") > 0) {
        config.stickiness = result["stickiness"].as<int>();
    }
#endif
    return config;
}

inline void print_pq_config(PriorityQueueConfig const& config) {
    std::clog << "factor: " << config.c << ", stickiness: " << config.stickiness;
}

template <typename PriorityQueue>
inline PriorityQueue create_pq(int num_threads, std::size_t initial_capacity, PriorityQueueConfig const& config) {
    return PriorityQueue(num_threads, initial_capacity, config);
}
#else
#ifdef PQ_CAPQ
#include "wrapper/capq.hpp"
using DefaultMinPriorityQueue = wrapper::Capq<true, true, true>;
static constexpr auto pq_name = "CA-PQ";
#elif defined PQ_KLSM
#ifndef KLSM_K
#error KLSM_K not defined
#endif
#include "wrapper/klsm.hpp"
#define PQ_HAS_GENERIC
template <typename Key, typename T>
using GenericMinPriorityQueue = wrapper::Klsm<Key, T, KLSM_K>;
using DefaultMinPriorityQueue = GenericMinPriorityQueue<unsigned long, unsigned long>;
static constexpr auto pq_name = "KLSM (k=" TOSTRING(KLSM_K) ")";
#elif defined PQ_LINDEN
#include "wrapper/linden.hpp"
using DefaultMinPriorityQueue = wrapper::Linden;
static constexpr auto pq_name = "Linden";
#elif defined PQ_SPRAYLIST
#include "wrapper/spraylist.hpp"
using DefaultMinPriorityQueue = wrapper::Spraylist;
static constexpr auto pq_name = "Spraylist";
#elif defined PQ_TBB_PQ
#include "wrapper/tbb_priority_queue.hpp"
#define PQ_HAS_GENERIC
#define PQ_HAS_MAX
#define PQ_HAS_MAX_GENERIC
template <typename Key, typename T>
using GenericMinPriorityQueue = wrapper::TBBPriorityQueue<Key, T, std::greater<>>;
template <typename Key, typename T>
using GenericMaxPriorityQueue = wrapper::TBBPriorityQueue<Key, T, std::less<>>;
using DefaultMinPriorityQueue = GenericMinPriorityQueue<unsigned long, unsigned long>;
using DefaultMaxPriorityQueue = GenericMaxPriorityQueue<unsigned long, unsigned long>;
static constexpr auto pq_name = "TBB";
#else
#error No valid PQ specified
#endif
#include "cxxopts.hpp"

inline void add_pq_options(cxxopts::Options&) {
}

struct PriorityQueueConfig {};

inline PriorityQueueConfig get_pq_options(cxxopts::ParseResult const&) {
    return PriorityQueueConfig{};
}

template <typename PriorityQueue>
inline PriorityQueue create_pq(int num_threads, std::size_t initial_capacity, PriorityQueueConfig const&) {
    return PriorityQueue(num_threads, initial_capacity);
}

inline void print_pq_config(PriorityQueueConfig const&) {
}
#endif

#undef STRINGIFY
#undef TOSTRING
