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
#error L1_CACHE_LINESIZE needs to be defined
#endif

namespace thread_coordination {

namespace detail {
template <class CharT, class Traits>
class ScopedOutput {
    using out_type = std::basic_ostream<CharT, Traits>;

   private:
    std::unique_lock<std::mutex> l_;
    out_type& out_;

   public:
    explicit ScopedOutput(std::mutex& mutex, out_type& out, unsigned int id)
        : out_{out} {
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
    utils::barrier barrier;
    alignas(L1_CACHE_LINESIZE) std::atomic_size_t index{0};
    std::chrono::steady_clock::time_point timestamp;
    std::mutex write_mutex;

    SharedData(unsigned int num_threads) : barrier{num_threads} {}
};

}  // namespace detail

class Context {
    detail::SharedData& shared_data_;
    unsigned int num_threads_;
    unsigned int id_;

   public:
    Context(detail::SharedData& sd, unsigned int nt, unsigned int id)
        : shared_data_{sd}, num_threads_{nt}, id_{id} {}

    unsigned int get_num_threads() const noexcept { return num_threads_; };

    unsigned int get_id() const noexcept { return id_; }

    template <class CharT, class Traits>
    detail::ScopedOutput<CharT, Traits> write(
        std::basic_ostream<CharT, Traits>& out) const noexcept {
        return detail::ScopedOutput<CharT, Traits>(shared_data_.write_mutex,
                                                   out, id_);
    }

    void synchronize() { shared_data_.barrier.wait(); }

    template <typename CompletionFunc>
    void synchronize(CompletionFunc&& f) {
        shared_data_.barrier.wait(std::forward<CompletionFunc>(f));
    }

    template <typename Work, typename... Args>
    void execute_synchronized(Work work, Args&&... args) {
        shared_data_.barrier.wait();
        work(std::forward<Args>(args)...);
        shared_data_.barrier.wait();
    }

    template <typename Work, typename... Args>
    void execute_synchronized_timed(
        std::chrono::steady_clock::duration& duration, Work work,
        Args&&... args) {
        shared_data_.barrier.wait([this] {
            shared_data_.timestamp = std::chrono::steady_clock::now();
        });
        work(std::forward<Args>(args)...);
        shared_data_.barrier.wait([this, &duration] {
            duration =
                std::chrono::steady_clock::now() - shared_data_.timestamp;
        });
    }

    template <typename Iter, typename Work>
    void execute_synchronized_blockwise(Iter begin, Iter end, Work work) {
        static constexpr std::size_t block_size = 1 << 8;

        auto n = static_cast<std::size_t>(end - begin);
        shared_data_.barrier.wait();
        while (true) {
            std::size_t block_begin = shared_data_.index.fetch_add(
                block_size, std::memory_order_relaxed);
            if (block_begin >= n) {
                break;
            }
            std::size_t block_end = std::min(block_begin + block_size, n);
            work(begin + block_begin, begin + block_end);
        }
        shared_data_.barrier.wait(
            [this] { shared_data_.index.store(0, std::memory_order_relaxed); });
    }

    template <typename Iter, typename Work>
    void execute_synchronized_blockwise_timed(
        std::chrono::steady_clock::duration& duration, Iter begin, Iter end,
        Work work) {
        static constexpr std::size_t block_size = 1 << 8;

        auto n = static_cast<std::size_t>(end - begin);
        shared_data_.barrier.wait([this] {
            shared_data_.timestamp = std::chrono::steady_clock::now();
        });
        while (true) {
            std::size_t block_begin = shared_data_.index.fetch_add(
                block_size, std::memory_order_relaxed);
            if (block_begin >= n) {
                break;
            }
            std::size_t block_end = std::min(block_begin + block_size, n);
            work(begin + block_begin, begin + block_end);
        }
        shared_data_.barrier.wait([this, &duration] {
            duration =
                std::chrono::steady_clock::now() - shared_data_.timestamp;
            shared_data_.index.store(0, std::memory_order_relaxed);
        });
    }
};

struct TaskHandle {
    detail::SharedData shared_data;
    std::vector<threading::pthread> threads;
    unsigned int num_threads;
    explicit TaskHandle(unsigned int n)
        : shared_data(n), num_threads{n} {}

    template <typename Task, typename... Args>
    void run_task(Args const&... args) {
        threads.clear();
        for (unsigned int i = 0; i < num_threads; ++i) {
            Context ctx{shared_data, num_threads, i};
            threads.emplace_back(Task::get_config(i), Task::run, ctx, args...);
        }
    }
    void join() {
        for (auto& t : threads) {
            t.join();
        }
    }
};

template <typename Task, typename... Args>
void run_task(unsigned int num_threads, Args const&... args) {
    auto shared_data = detail::SharedData(num_threads);
    std::vector<threading::pthread> threads;
    for (unsigned int i = 0; i < num_threads; ++i) {
        Context ctx{shared_data, num_threads, i};
        threads.emplace_back(Task::get_config(i), Task::run, ctx, args...);
    }
    for (auto& t : threads) {
        t.join();
    }
}

}  // namespace thread_coordination

#endif
