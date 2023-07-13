#pragma once

#include "barrier.hpp"
#include "threading.hpp"

#include <pthread.h>

#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <ostream>
#include <utility>
#include <vector>

namespace thread_coordination {

using clock_type = std::chrono::steady_clock;
struct time_result_type {
    clock_type::time_point start;
    clock_type::time_point end;
};

namespace detail {
template <class CharT, class Traits>
class ScopedOutput {
    using out_type = std::basic_ostream<CharT, Traits>;

   private:
    std::unique_lock<std::mutex> l_;
    out_type& out_;

   public:
    explicit ScopedOutput(std::mutex& mutex, out_type& out, int id) : out_{out} {
        l_ = std::unique_lock{mutex};
        out_ << "[Thread " << id << "] ";
    }

    template <typename T>
    friend out_type& operator<<(ScopedOutput const& so, T const& t) {
        so.out_ << t;
        return so.out_;
    }
};

struct SharedData {
    int num_threads;
    utils::Barrier barrier;
    std::mutex write_mutex;
    alignas(64) std::atomic_llong counter{0};

    explicit SharedData(int n) : num_threads{n}, barrier{n} {
    }
};

}  // namespace detail

class Context {
    template <typename Affinity>
    friend class TaskHandle;

    detail::SharedData& shared_data_;
    int id_;

    Context(detail::SharedData& sd, int id) : shared_data_{sd}, id_{id} {
    }

   public:
    [[nodiscard]] int get_num_threads() const noexcept {
        return shared_data_.num_threads;
    };

    [[nodiscard]] int get_id() const noexcept {
        return id_;
    }

    template <class CharT, class Traits>
    detail::ScopedOutput<CharT, Traits> write(std::basic_ostream<CharT, Traits>& out) const {
        return detail::ScopedOutput<CharT, Traits>(shared_data_.write_mutex, out, id_);
    }

    void synchronize() const {
        shared_data_.barrier.wait();
    }

    template <typename CompletionFunc>
    void synchronize(CompletionFunc&& f) const {
        shared_data_.barrier.wait(std::forward<CompletionFunc>(f));
    }

    template <typename Work, typename... Args>
    time_result_type execute_synchronized(Work work, Args&&... args) const {
        shared_data_.barrier.wait();
        auto t_start = clock_type::now();
        work(std::forward<Args>(args)...);
        auto t_end = clock_type::now();
        shared_data_.barrier.wait();
        return {t_start, t_end};
    }

    template <typename Work, typename... Args>
    time_result_type execute_synchronized_blockwise(long long n, Work work, Args&&... args) const {
        static constexpr auto block_size = static_cast<long long>(1) << 12;

        shared_data_.barrier.wait();
        auto t_start = clock_type::now();
        while (true) {
            auto begin = shared_data_.counter.fetch_add(block_size, std::memory_order_relaxed);
            if (begin >= n) {
                break;
            }
            auto end = std::min(n, begin + block_size);
            work(begin, end, args...);  // no perfect forwarding here
        }
        auto t_end = clock_type::now();
        shared_data_.barrier.wait([this] { shared_data_.counter.store(0, std::memory_order_relaxed); });
        return {t_start, t_end};
    }
};

namespace affinity {

struct individual_cores {
    std::size_t stride = 1;
    std::size_t offset = 0;

    threading::thread_config operator()(int id) const {
        threading::thread_config cfg;
        cfg.cpu_set.set(offset + static_cast<std::size_t>(id) * stride);
        return cfg;
    }
};

struct same_core {
    std::size_t core = 1;

    threading::thread_config operator()(int /*unused*/) const {
        threading::thread_config cfg;
        cfg.cpu_set.set(core);
        return cfg;
    }
};

struct NUMA {
    int cores_per_node = 4;
    int num_nodes = 16;

    threading::thread_config operator()(int id) const {
        threading::thread_config cfg;
        int cpu = 0;
        if (id % (2 * cores_per_node) < cores_per_node) {
            cpu = (id / cores_per_node) * (cores_per_node / 2) + (id % cores_per_node);
        } else {
            id -= cores_per_node;
            cpu = num_nodes * cores_per_node + (id / cores_per_node) * cores_per_node / 2 + (id % cores_per_node);
        }
        cfg.cpu_set.set(static_cast<std::size_t>(cpu));
        return cfg;
    }
};

}  // namespace affinity

template <typename Affinity = affinity::individual_cores>
class TaskHandle {
    std::vector<threading::pthread> threads_;
    detail::SharedData shared_data_;

   public:
    template <typename Task, typename... Args>
    explicit TaskHandle(Affinity affinity, int num_threads, Task task, Args... args)
        : threads_(static_cast<std::size_t>(num_threads)), shared_data_(num_threads) {
        for (std::size_t i = 0; i < threads_.size(); ++i) {
            Context ctx{shared_data_, static_cast<int>(i)};
            threads_[i] = threading::pthread(affinity(static_cast<int>(i)), task, ctx, args...);
        }
    }

    template <typename Task, typename... Args>
    explicit TaskHandle(int num_threads, Task task, Args... args) : TaskHandle(Affinity{}, num_threads, task, args...) {
    }

    TaskHandle(TaskHandle const&) = delete;
    TaskHandle(TaskHandle&&) noexcept = default;
    TaskHandle& operator=(TaskHandle const&) = delete;
    TaskHandle& operator=(TaskHandle&&) noexcept = default;
    ~TaskHandle() = default;

    void wait() {
        for (auto& t : threads_) {
            t.join();
        }
    }
};

}  // namespace thread_coordination
