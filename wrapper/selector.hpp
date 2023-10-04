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

#include "wrapper/priority.hpp"

#if defined PQ_MQ
#include "wrapper/multiqueue.hpp"
template <typename Key, typename T, Priority P>
using PriorityQueue = wrapper::MultiQueue<Key, T, P>;
#elif defined PQ_CAPQ
#include "wrapper/capq.hpp"
template <typename Key, typename T, Priority P>
using PriorityQueue = wrapper::CAPQ<Key, T, P>;
#elif defined PQ_KLSM
#include "wrapper/klsm.hpp"
template <typename Key, typename T, Priority P>
using PriorityQueue = wrapper::KLSM<Key, T, P>;
#elif defined PQ_LINDEN
#include "wrapper/linden.hpp"
template <typename Key, typename T, Priority P>
using PriorityQueue = wrapper::Linden<Key, T, P>;
#elif defined PQ_SPRAYLIST
#include "wrapper/spraylist.hpp"
template <typename Key, typename T, Priority P>
using PriorityQueue = wrapper::Spraylist<Key, T, P>;
#elif defined PQ_TBB_PQ
#include "wrapper/tbb_priority_queue.hpp"
template <typename Key, typename T, Priority P>
using PriorityQueue =
    wrapper::TBBPriorityQueue<Key, T, P>;
#elif defined PQ_TBB_FIFO
#include "wrapper/tbb_queue.hpp"
template <typename Key, typename T, Priority P>
using PriorityQueue = wrapper::TBBQueue<Key, T, P>;
#elif defined PQ_SMQ
#include "wrapper/smq.hpp"
template <typename Key, typename T, Priority P>
using PriorityQueue = wrapper::StealingMQ<Key, T, P>;
#else
#error No valid PQ specified
#endif
