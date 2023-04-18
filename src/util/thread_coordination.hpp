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
using timepoint_type = clock_type::time_point;

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
    using clock_type = std::chrono::steady_clock;

    int num_threads;
    utils::Barrier barrier;
    std::mutex write_mutex;
    alignas(64) std::atomic_size_t index{0};

    explicit SharedData(int n) : num_threads{n}, barrier{n} {
    }
};

}  // namespace detail

class Context {
    template <typename Affinity>
    friend class TaskHandle;

    using duration_type = detail::SharedData::clock_type::duration;

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
    std::pair<timepoint_type, timepoint_type> execute_synchronized(Work work, Args&&... args) const {
        shared_data_.barrier.wait();
        auto t1 = std::chrono::steady_clock::now();
        work(std::forward<Args>(args)...);
        auto t2 = std::chrono::steady_clock::now();
        shared_data_.barrier.wait();
        return {t1, t2};
    }

    template <typename Iter, typename Work, typename... Args>
    std::pair<timepoint_type, timepoint_type> execute_synchronized_blockwise(Iter begin, Iter end, Work work, Args&&... args) const {
        using difference_type = typename std::iterator_traits<Iter>::difference_type;
        static constexpr difference_type block_size = 1 << 12;

        difference_type n = end - begin;
        shared_data_.barrier.wait();
        auto t1 = std::chrono::steady_clock::now();
        while (true) {
            auto block_begin =
                static_cast<difference_type>(shared_data_.index.fetch_add(block_size, std::memory_order_relaxed));
            if (block_begin >= n) {
                break;
            }
            auto block_end = std::min(block_begin + block_size, n);
            work(begin + block_begin, begin + block_end, args...); // no perfect forwarding here
        }
        auto t2 = std::chrono::steady_clock::now();
        shared_data_.barrier.wait([this] { shared_data_.index.store(0, std::memory_order_relaxed); });
        return {t1, t2};
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
