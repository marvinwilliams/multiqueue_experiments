#pragma once

#ifdef __SSE2__
#include <emmintrin.h>
#define PAUSE _mm_pause()
#else
#define PAUSE void(0)
#endif
#include <atomic>

namespace termination_detection {

class TerminationDetection {
    int num_threads_;
    std::atomic_int idle_count_{0};
    std::atomic_int no_work_count_{0};

    bool should_terminate() {
        idle_count_.fetch_add(1, std::memory_order_relaxed);
        while (no_work_count_.load(std::memory_order_relaxed) >= num_threads_) {
            if (idle_count_.load(std::memory_order_relaxed) >= num_threads_) {
                return true;
            }
            PAUSE;
        }
        idle_count_.fetch_sub(1, std::memory_order_relaxed);
        return false;
    }

   public:
    explicit TerminationDetection(int num_threads) : num_threads_{num_threads} {
    }

    void reset() {
        idle_count_.store(0, std::memory_order_relaxed);
        no_work_count_.store(0, std::memory_order_relaxed);
    }

    template <typename F>
    bool repeat(F&& f) {
        if (f()) {
            return true;
        }
        no_work_count_.fetch_add(1, std::memory_order_relaxed);
        while (!f()) {
            if (no_work_count_.load(std::memory_order_relaxed) >= num_threads_) {
                if (should_terminate()) {
                    return false;
                }
            }
        }
        no_work_count_.fetch_sub(1, std::memory_order_relaxed);
        return true;
    }
};

}  // namespace termination_detection
