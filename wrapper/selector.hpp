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

#include <utility>

#if defined PQ_MQ
#include "wrapper/multiqueue.hpp"
using namespace wrapper::multiqueue;
#elif defined PQ_CAPQ
#include "wrapper/capq.hpp"
using namespace wrapper::capq;
#elif defined PQ_KLSM
#include "wrapper/klsm.hpp"
using namespace wrapper::klsm;
#elif defined PQ_LINDEN
#include "wrapper/linden.hpp"
using namespace wrapper::linden;
#elif defined PQ_SPRAYLIST
#include "wrapper/spraylist.hpp"
using namespace wrapper::spraylist;
#elif defined PQ_TBB_PQ
#include "wrapper/tbb_priority_queue.hpp"
using namespace wrapper::TBBPriorityQueue;
#elif defined PQ_TBB_FIFO
#include "wrapper/tbb_queue.hpp"
using namespace wrapper::TBBQueue;
#elif defined PQ_SMQ
#include <functional>
#include "wrapper/smq.hpp"
using namespace wrapper::smq;
#else
#error No valid PQ specified
#endif
