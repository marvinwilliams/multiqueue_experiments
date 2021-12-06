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

#if defined PQ_MQ_RANDOM || defined PQ_MQ_STICKY
#define PQ_IS_MQ
#else
#define PQ_IS_WRAPPER
#endif

#if defined PQ_IS_MQ
#include "multiqueue/default_configuration.hpp"
#include "multiqueue/multiqueue.hpp"
#include "multiqueue/selection_strategy/random.hpp"
#include "multiqueue/selection_strategy/sticky.hpp"
#elif defined PQ_CAPQ || defined PQ_CAPQ1 || defined PQ_CAPQ2 || \
    defined PQ_CAPQ3 || defined PQ_CAPQ4
#include "capq.hpp"
#elif defined PQ_KLSM || defined PQ_KLSM256 || defined PQ_KLSM1024 || \
    defined PQ_KLSM4096
#include "klsm.hpp"
#elif defined PQ_LINDEN
#include "linden.hpp"
#elif defined PQ_SPRAYLIST
#include "spraylist.hpp"
#elif defined PQ_FIFO

#elif defined PQ_LOCKING

#elif defined PQ_TBB

#else
#error No valid PQ specified
#endif

#include <cstdint>
#include <functional>
#include <sstream>
#include <string>
#include <type_traits>

namespace util {

#if defined PQ_IS_MQ

struct Config : multiqueue::DefaultConfiguration {
#ifdef PQ_MQ_RANDOM
  using selection_strategy_t = multiqueue::selection_strategy::random;
#elif PQ_MQ_STICKY
  using selection_strategy_t = multiqueue::selection_strategy::sticky;
#endif
#ifdef MQ_NO_BUFFERING
  static constexpr bool UseBuffers = false;
#endif
#ifdef MQ_DBUF_SIZE
  static constexpr std::size_t DeletionBufferSize = MQ_DBUF_SIZE;
#endif
#ifdef MQ_IBUF_SIZE
  static constexpr std::size_t InsertionBufferSize = MQ_IBUF_SIZE;
#endif
#ifdef MQ_HEAP_DEGREE
  static constexpr unsigned int HeapDegree = MQ_HEAP_DEGREE;
#endif
#ifdef MQ_IMPLICIT_LOCK
  static constexpr bool ImplicitLock = true;
#endif
};
#endif

template <typename KeyType, typename ValueType>
struct PriorityQueueFactory {
#if defined PQ_IS_MQ
  using type = typename multiqueue::MultiqueueFactory<
      KeyType, ValueType, std::less<>, Config, std::allocator<KeyType>>::type;
#elif defined PQ_CAPQ || defined PQ_CAPQ1 || defined PQ_CAPQ2 || \
    defined PQ_CAPQ3 || defined PQ_CAPQ4
  // not available with generic types
#elif defined PQ_KLSM256
  using type = wrapper::Klsm<KeyType, ValueType, 256>;
#elif defined PQ_KLSM1024
  using type = wrapper::Klsm<KeyType, ValueType, 1024>;
#elif defined PQ_KLSM4096
  using type = wrapper::Klsm<KeyType, ValueType, 4096>;
#elif defined PQ_LINDEN
  // not available with generic types
#elif defined PQ_SPRAYLIST
  // not available with generic types
#endif
};

template <>
struct PriorityQueueFactory<unsigned long, unsigned long> {
  using KeyType = unsigned long;
  using ValueType = unsigned long;
#if defined PQ_IS_MQ
  using type = typename multiqueue::MultiqueueFactory<
      KeyType, ValueType, std::less<>, Config, std::allocator<KeyType>>::type;
#elif defined PQ_CAPQ || defined PQ_CAPQ1
  using type = wrapper::Capq<true, true, true>;
#elif defined PQ_CAPQ2
  using type = wrapper::Capq<true, false, true>;
#elif defined PQ_CAPQ3
  using type = wrapper::Capq<false, true, true>;
#elif defined PQ_CAPQ4
  using type = wrapper::Capq<false, false, true>;
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
#endif
};

}  // namespace util

#endif  //! UTILS_PRIORITY_QUEUE_FACTORY_HPP_INCLUDED
