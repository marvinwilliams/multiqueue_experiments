#pragma once

#include "util/benchmark.hpp"
#include "util/thread_coordination.hpp"

#include <atomic>
#include <optional>

namespace parallel_search {

class Termination {
    int num_threads_;
    std::atomic_int idle_count_{0};
    std::atomic_int no_work_count_{0};
    std::atomic_llong missing_nodes_{0};

    bool should_terminate() {
        idle_count_.fetch_add(1, std::memory_order_relaxed);
        while (no_work_count_.load(std::memory_order_relaxed) >= num_threads_) {
            if (idle_count_.load(std::memory_order_relaxed) >= num_threads_) {
                return true;
            }
            thread_coordination::cpu_relax();
        }
        idle_count_.fetch_sub(1, std::memory_order_relaxed);
        return false;
    }

    template <typename F>
    bool repeat(F&& f) {
        if (f()) {
            return true;
        }
        no_work_count_.fetch_add(1, std::memory_order_relaxed);
        while (!f()) {
            if (no_work_count_.load(std::memory_order_relaxed) >= num_threads_ && should_terminate()) {
                return false;
            }
        }
        no_work_count_.fetch_sub(1, std::memory_order_relaxed);
        return true;
    }

   public:
    explicit Termination(int num_threads) noexcept : num_threads_{num_threads} {
    }

    template <typename Handle, typename Process, typename NodeCount>
    benchmark::Interval run(thread_coordination::Context& ctx, Handle& handle, Process&& process,
                            NodeCount&& node_count) {
        benchmark::Interval interval;
        ctx.synchronize();
        ctx.spin_synchronize();
        interval.start = benchmark::clock_type::now();
        while (true) {
            std::optional<typename Handle::value_type> node;
            while (repeat([&]() {
                node = handle.try_pop();
                return node.has_value();
            })) {
                process(*node);
            }
            missing_nodes_.fetch_add(node_count(), std::memory_order_relaxed);
            ctx.spin_synchronize();
            if (missing_nodes_.load(std::memory_order_relaxed) == 0) {
                interval.end = benchmark::clock_type::now();
                return interval;
            }
            ctx.spin_synchronize();
            if (ctx.id() == 0) {
                missing_nodes_.store(0, std::memory_order_relaxed);
                idle_count_.store(0, std::memory_order_relaxed);
                no_work_count_.store(0, std::memory_order_relaxed);
            }
            ctx.spin_synchronize();
        }
    }
};

}  // namespace parallel_search
