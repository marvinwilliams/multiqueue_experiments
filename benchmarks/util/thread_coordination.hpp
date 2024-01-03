#pragma once

#include "affinity.hpp"
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

class Context {
    friend class Runner;

    int id_;
    Barrier* barrier_;
    std::mutex* write_mutex_{};

    Context(int id, Barrier& barrier, std::mutex& m) : id_{id}, barrier_{&barrier}, write_mutex_{&m} {
    }

    class GuardedWriter {
        friend Context;
        std::unique_lock<std::mutex> lock_;
        std::ostream* out_;

        explicit GuardedWriter(std::mutex& mutex, std::ostream& out, int id) : lock_{mutex}, out_{&out} {
            *out_ << "[Thread " << id << "] ";
        }

       public:
        template <typename T>
        std::ostream& operator<<(T&& t) {
            *out_ << std::forward<T>(t);
            return *out_;
        }
    };

   public:
    Context(Context const&) = delete;
    Context(Context&&) noexcept = default;

    Context& operator=(Context const&) = delete;
    Context& operator=(Context&&) noexcept = default;

    ~Context() = default;

    [[nodiscard]] int id() const noexcept {
        return id_;
    }

    GuardedWriter write(std::ostream& out) {
        return GuardedWriter(*write_mutex_, out, id_);
    }

    void synchronize() const {
        barrier_->wait();
    }

    template <typename CompletionFunc>
    void synchronize(CompletionFunc&& f) const {
        barrier_->wait(std::forward<CompletionFunc>(f));
    }
};

class Runner {
    std::vector<threading::pthread> threads_;
    Barrier barrier_;
    std::mutex write_mutex_;

   public:
    template <typename Affinity, typename Task, typename... Args>
    explicit Runner(Affinity affinity, int num_threads, Task task, Args... args)
        : threads_(static_cast<std::size_t>(num_threads)), barrier_(num_threads) {
        for (int i = 0; i < num_threads; ++i) {
            Context context{i, barrier_, write_mutex_};
            threads_[i] = threading::pthread(affinity(i), task, context, args...);
        }
    }

    template <typename Task, typename... Args>
    explicit Runner(int num_threads, Task task, Args&&... args)
        : Runner(affinity::individual_cores{}, num_threads, task, std::forward<Args>(args)...) {
    }

    void wait() {
        for (auto& t : threads_) {
            t.join();
        }
    }
};

}  // namespace thread_coordination
