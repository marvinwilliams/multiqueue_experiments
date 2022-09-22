#ifndef THREAD_COORDINATION_HPP_7F7173D6
#define THREAD_COORDINATION_HPP_7F7173D6

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

#ifndef L1_CACHE_LINESIZE
#define L1_CACHE_LINESIZE 64
#endif

namespace thread_coordination {

using duration_type = std::chrono::steady_clock::duration;

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
    utils::barrier barrier;
    alignas(L1_CACHE_LINESIZE) std::atomic_size_t index{0};
    clock_type::time_point timestamp;
    std::mutex mutex;
    std::mutex write_mutex;

    explicit SharedData(int num_threads) : barrier{num_threads} {
    }
};

}  // namespace detail

class Context {
    friend class TaskHandle;
    using duration_type = detail::SharedData::clock_type::duration;

    detail::SharedData& shared_data_;
    int num_threads_;
    int id_;

    Context(detail::SharedData& sd, int n, int id) : shared_data_{sd}, num_threads_{n}, id_{id} {
    }

   public:
    [[nodiscard]] int get_num_threads() const noexcept {
        return num_threads_;
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
    void execute_synchronized(Work work, Args&&... args) const {
        shared_data_.barrier.wait();
        work(std::forward<Args>(args)...);
        shared_data_.barrier.wait();
    }

    template <typename Work, typename... Args>
    void execute_synchronized_timed(duration_type& duration, Work work, Args&&... args) const {
        shared_data_.barrier.wait([this] { shared_data_.timestamp = std::chrono::steady_clock::now(); });
        work(std::forward<Args>(args)...);
        shared_data_.barrier.wait(
            [this, &duration] { duration = std::chrono::steady_clock::now() - shared_data_.timestamp; });
    }

    template <typename Iter, typename Work>
    void execute_synchronized_blockwise(Iter begin, Iter end, Work work) const {
        static constexpr std::size_t block_size = 1 << 12;

        auto n = static_cast<std::size_t>(end - begin);
        shared_data_.barrier.wait();
        while (true) {
            std::size_t block_begin = shared_data_.index.fetch_add(block_size, std::memory_order_relaxed);
            if (block_begin >= n) {
                break;
            }
            std::size_t block_end = std::min(block_begin + block_size, n);
            work(begin + block_begin, begin + block_end);
        }
        shared_data_.barrier.wait([this] { shared_data_.index.store(0, std::memory_order_relaxed); });
    }

    template <typename Iter, typename Work>
    void execute_synchronized_blockwise_timed(duration_type& duration, Iter begin, Iter end, Work work) const {
        static constexpr std::size_t block_size = 1 << 12;

        auto n = static_cast<std::size_t>(end - begin);
        shared_data_.barrier.wait([this] { shared_data_.timestamp = std::chrono::steady_clock::now(); });
        while (true) {
            std::size_t block_begin = shared_data_.index.fetch_add(block_size, std::memory_order_relaxed);
            if (block_begin >= n) {
                break;
            }
            std::size_t block_end = std::min(block_begin + block_size, n);
            work(begin + block_begin, begin + block_end);
        }
        shared_data_.barrier.wait([this, &duration] {
            duration = std::chrono::steady_clock::now() - shared_data_.timestamp;
            shared_data_.index.store(0, std::memory_order_relaxed);
        });
    }

    template <typename Func>
    auto execute_exclusive(Func f) {
        std::scoped_lock l{shared_data_.mutex};
        return f();
    }
};

class TaskHandle {
    detail::SharedData shared_data;
    std::vector<threading::pthread> threads;

   public:
    template <typename Task, typename... Args>
    explicit TaskHandle(int num_threads, Args const&... args) : shared_data(num_threads) {
        for (int i = 0; i < num_threads; ++i) {
            Context ctx{shared_data, num_threads, i};
            threads.emplace_back(Task::get_config(i), Task::run, ctx, args...);
        }
    }

    void join() {
        threads.clear();
    }
};

}  // namespace thread_coordination

#endif
