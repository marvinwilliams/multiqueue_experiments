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
#include "multiqueue/multiqueue.hpp"
#include "multiqueue/queue_selection/global_permutation.hpp"
#include "multiqueue/queue_selection/random.hpp"
#include "multiqueue/queue_selection/stick_random.hpp"
#include "multiqueue/queue_selection/swap_assignment.hpp"

#include "cxxopts.hpp"

#include <iostream>
#ifdef USE_STD_PQ
#include <queue>
#endif

#define PQ_HAS_GENERIC
#define PQ_HAS_MAX
#define PQ_HAS_MAX_GENERIC
namespace detail {
template <typename Key, typename T, typename Compare>
using GenericPriorityQueue = multiqueue::MultiQueue<
    Key, std::pair<Key, T>, Compare
#ifdef USE_STD_PQ
    ,
    multiqueue::defaults::MultiQueueTraits, multiqueue::defaults::KeyOfValue<Key, std::pair<Key, T>>,
    multiqueue::BufferedPQ<
        std::priority_queue<std::pair<Key, T>, std::vector<std::pair<Key, T>>,
                            multiqueue::defaults::ValueCompare<
                                std::pair<Key, T>, multiqueue::defaults::KeyOfValue<Key, std::pair<Key, T>>, Compare>>>
#endif
    >;
}  // namespace detail

template <typename Key, typename T>
using GenericMinPriorityQueue = detail::GenericPriorityQueue<Key, T, std::greater<>>;
template <typename Key, typename T>
using GenericMaxPriorityQueue = detail::GenericPriorityQueue<Key, T, std::less<>>;
using DefaultMinPriorityQueue = GenericMinPriorityQueue<unsigned long, unsigned long>;
using DefaultMaxPriorityQueue = GenericMaxPriorityQueue<unsigned long, unsigned long>;
inline void add_pq_options(cxxopts::Options& options) {
    options.add_options()("c,factor", "The factor for queues", cxxopts::value<int>(), "NUMBER")
#ifndef PQ_MQ_NONE
        ("k,stickiness", "The stickiness period", cxxopts::value<int>(), "NUMBER");
#endif
}

struct PriorityQueueConfig {
    multiqueue::build_config::DefaultQueueSelectionPolicy::Config queue_selection_policy_config;
    int pqs_per_thread = 2;
};

inline PriorityQueueConfig get_pq_options(cxxopts::ParseResult const& result) {
    PriorityQueueConfig config;
    if (result.count("factor") > 0) {
        config.pqs_per_thread = result["factor"].as<int>();
    }
#ifndef PQ_MQ_NONE
    if (result.count("stickiness") > 0) {
        config.queue_selection_policy_config.stickiness = result["stickiness"].as<int>();
    }
#endif
    return config;
}

template <typename Policy>
struct QueueSelectionPolicyName {
    static constexpr auto name = "Unknown";
};

template <>
struct QueueSelectionPolicyName<multiqueue::queue_selection::Random> {
    static constexpr auto name = "Random";
};

template <>
struct QueueSelectionPolicyName<multiqueue::queue_selection::StickRandom> {
    static constexpr auto name = "StickRandom";
};

template <>
struct QueueSelectionPolicyName<multiqueue::queue_selection::SwapAssignment> {
    static constexpr auto name = "Swapping";
};

template <>
struct QueueSelectionPolicyName<multiqueue::queue_selection::GlobalPermutation> {
    static constexpr auto name = "GlobalPermutation";
};

inline void describe_pq(std::ostream& out, PriorityQueueConfig const& config) {
    out << "MultiQueue (" << QueueSelectionPolicyName<multiqueue::build_config::DefaultQueueSelectionPolicy>::name
        << ", heap arity: " << multiqueue::build_config::DefaultHeapArity
        << ", buffer sizes (ins/del): " << multiqueue::build_config::DefaultInsertionBuffersize << '/'
        << multiqueue::build_config::DefaultDeletionBuffersize << ", pqs per thread: " << config.pqs_per_thread
#ifndef PQ_MQ_NONE
        << ", stickiness = " << config.queue_selection_policy_config.stickiness
#endif
        << ')';
}

template <typename PriorityQueue>
inline PriorityQueue create_pq(int num_threads, std::size_t initial_capacity, PriorityQueueConfig const& config) {
    return PriorityQueue(num_threads * config.pqs_per_thread, initial_capacity, config.queue_selection_policy_config);
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
static constexpr auto pq_name = "KLSM (k: " TOSTRING(KLSM_K) ")";
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

inline void describe_pq(std::ostream& out, PriorityQueueConfig const&) {
    out << pq_name;
}

template <typename PriorityQueue>
inline PriorityQueue create_pq(int num_threads, std::size_t initial_capacity, PriorityQueueConfig const&) {
    return PriorityQueue(num_threads, initial_capacity);
}

#endif

#undef STRINGIFY
#undef TOSTRING
