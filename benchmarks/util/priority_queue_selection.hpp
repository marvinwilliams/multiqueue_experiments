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

#include <functional>
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
