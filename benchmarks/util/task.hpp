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

namespace task {

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

    explicit SharedData(int n) : num_threads{n}, barrier{n} {
    }
};

}  // namespace detail

class Control {
    template <typename Affinity>
    friend class Runner;

    detail::SharedData& shared_data_;
    const int id_;

    Control(detail::SharedData& sd, int id) : shared_data_{sd}, id_{id} {
    }

   public:
    [[nodiscard]] int num_threads() const noexcept {
        return shared_data_.num_threads;
    };

    [[nodiscard]] int id() const noexcept {
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

    template <typename Func>
    void once(Func f) const {
        if (id_ == 0) {
            f();
        }
    }
};

template <typename Affinity = affinity::individual_cores>
class Runner {
    detail::SharedData shared_data_;
    std::vector<threading::pthread> threads_;

   public:
    template <typename Task, typename... Args>
    explicit Runner(Affinity affinity, int num_threads, Task task, Args... args)
        : shared_data_(num_threads), threads_(static_cast<std::size_t>(num_threads)) {
        for (std::size_t i = 0; i < threads_.size(); ++i) {
            Control control{shared_data_, static_cast<int>(i)};
            threads_[i] = threading::pthread(affinity(static_cast<int>(i)), task, control, args...);
        }
    }

    template <typename Task, typename... Args>
    explicit Runner(int num_threads, Task task, Args... args) : Runner(Affinity{}, num_threads, task, args...) {
    }

    Runner(Runner const&) = delete;
    Runner(Runner&&) noexcept = default;
    Runner& operator=(Runner const&) = delete;
    Runner& operator=(Runner&&) noexcept = default;
    ~Runner() = default;  // must wait before destruction

    void wait() {
        for (auto& t : threads_) {
            t.join();
        }
    }
};

}  // namespace task
