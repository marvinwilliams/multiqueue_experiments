/**
******************************************************************************
* @file:   selector.hpp
*
* @author: Marvin Williams
* @date:   2021/02/22 13:23
* @brief:
*******************************************************************************
**/
#pragma once

#if defined PQ_MQ
#include "wrapper/multiqueue.hpp"
#define PQ wrapper::multiqueue::MultiQueue
#elif defined PQ_CAPQ
#include "wrapper/capq.hpp"
#define PQ wrapper::capq::CAPQ
#elif defined PQ_KLSM
#include "wrapper/klsm.hpp"
#define PQ wrapper::klsm::KLsm
#elif defined PQ_LINDEN
#include "wrapper/linden.hpp"
#define PQ wrapper::linden::Linden
#elif defined PQ_SPRAYLIST
#include "wrapper/spraylist.hpp"
#define PQ wrapper::spraylist::Spraylist
#elif defined PQ_TBB_PQ
#include "wrapper/tbb_priority_queue.hpp"
#define PQ wrapper::tbb_pq::TBBPriorityQueue
#elif defined PQ_TBB_FIFO
#include "wrapper/tbb_queue.hpp"
#define PQ wrapper::tbb_fifo::TBBFIFO
#elif defined PQ_SMQ
#include "wrapper/smq.hpp"
#define PQ wrapper::stealing_mq::StealingMQ
#elif defined PQ_LOCKED_PQ
#include "wrapper/locked_pq.hpp"
#define PQ wrapper::locked_pq::LockedPQ
#elif defined PQ_MQ_PQ
#include "wrapper/mq_pq.hpp"
#define PQ wrapper::mq_pq::MultiQueuePQ
#else
#error No valid PQ specified
#endif
