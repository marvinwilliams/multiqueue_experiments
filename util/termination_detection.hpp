#pragma once

#ifdef __SSE2__
#include <emmintrin.h>
#define PAUSE _mm_pause()
#else
#define PAUSE void(0)
#endif
#include <atomic>
#include <condition_variable>
#include <mutex>

namespace termination_detection {

struct Data {
    std::atomic_int idle_count{0};
    std::atomic_int no_work_count{0};
};

namespace detail {

inline bool wait_to_terminate(int num_threads, Data& data) {
    auto idle_count = data.idle_count.fetch_add(1, std::memory_order_relaxed) + 1;
    while (idle_count < num_threads) {
        if (data.no_work_count.load(std::memory_order_relaxed) < num_threads) {
            data.idle_count.fetch_sub(1, std::memory_order_relaxed);
            return false;
        }
        PAUSE;
        idle_count = data.idle_count.load(std::memory_order_relaxed);
    }
    return true;
}

}  // namespace detail

template <typename F>
bool try_do(int num_threads, Data& data, F f) {
    for (int i = 0; i < 100; ++i) {
        if (f()) {
            return true;
        }
    }
    auto num_no_work = data.no_work_count.fetch_add(1, std::memory_order_relaxed) + 1;
    while (!f()) {
        if (num_no_work >= num_threads && detail::wait_to_terminate(num_threads, data)) {
            return false;
        }
        num_no_work = data.no_work_count.load(std::memory_order_relaxed);
    }
    data.no_work_count.fetch_sub(1, std::memory_order_relaxed);
    return true;
}

}  // namespace termination_detection
