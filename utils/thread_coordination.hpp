#ifndef THREAD_COORDINATION_HPP_7F7173D6
#define THREAD_COORDINATION_HPP_7F7173D6

#include "barrier.hpp"
#include "threading.hpp"

#include <pthread.h>

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <unordered_map>
#include <utility>
#include <vector>

namespace thread_coordination {

class ThreadCoordinator;

class Context {
  friend ThreadCoordinator;
  ThreadCoordinator& coordinator_;
  unsigned int id_;

 private:
  Context(ThreadCoordinator& coordinator, unsigned int id)
      : coordinator_{coordinator}, id_{id} {}

 public:
  bool is_main() const noexcept { return id_ == 0; }

  unsigned int get_id() const noexcept { return id_; }

  unsigned int get_num_threads() const noexcept;

  void synchronize();

  template <typename CompletionFunc>
  void synchronize(CompletionFunc&& f);

  template <typename Work, typename... Args>
  void execute_synchronized(Work work, Args&&... args);

  template <typename Work, typename... Args>
  void execute_synchronized_timed(std::string name, Work work, Args&&... args);

  template <typename Iter, typename Work>
  void execute_synchronized_blockwise(Iter begin, Iter end, Work work);

  template <typename Iter, typename Work>
  void execute_synchronized_blockwise_timed(std::string name, Iter begin,
                                            Iter end, Work work);
};

class ThreadCoordinator {
  friend Context;
  static constexpr std::ptrdiff_t block_size = 4096;

  unsigned int num_threads_;
  utils::barrier barrier_;
  std::vector<threading::pthread> threads_;
  alignas(L1_CACHE_LINESIZE) std::atomic_size_t index_;
  std::chrono::steady_clock::time_point start_;
  std::unordered_map<std::string, std::chrono::steady_clock::duration>
      durations_;

  void reset_index() { index_.store(0, std::memory_order_relaxed); }
  std::size_t get_next_index() {
    return index_.fetch_add(block_size, std::memory_order_relaxed);
  }

 public:
  explicit ThreadCoordinator(unsigned int num_threads)
      : num_threads_{num_threads}, barrier_{num_threads} {}

  unsigned int get_num_threads() const noexcept { return num_threads_; }

  template <typename Task, typename... Args>
  void run_task(Args const&... args) {
    for (unsigned int i = 0; i < num_threads_; ++i) {
      Context ctx{*this, i};
      threading::thread_config config = Task::get_config(ctx);
      threads_.emplace_back(config, Task::run, ctx, args...);
    }
  }

  void join() {
    for (auto& t : threads_) {
      if (t.joinable()) {
        t.join();
      }
    }
  }

  std::optional<std::chrono::steady_clock::duration> get_duration(
      std::string const& s) {
    if (auto it = durations_.find(s); it != durations_.end()) {
      return it->second;
    }
    return {};
  }
};

inline unsigned int Context::get_num_threads() const noexcept {
  return coordinator_.get_num_threads();
}

inline void Context::synchronize() { coordinator_.barrier_.wait(); }

template <typename CompletionFunc>
inline void Context::synchronize(CompletionFunc&& f) {
  coordinator_.barrier_.wait(std::forward<CompletionFunc>(f));
}

template <typename Work, typename... Args>
inline void Context::execute_synchronized(Work work, Args&&... args) {
  coordinator_.barrier_.wait();
  work(std::forward<Args>(args)...);
  coordinator_.barrier_.wait();
}

template <typename Work, typename... Args>
inline void Context::execute_synchronized_timed(std::string name, Work work,
                                                Args&&... args) {
  coordinator_.barrier_.wait(
      [this] { coordinator_.start_ = std::chrono::steady_clock::now(); });
  work(std::forward<Args>(args)...);
  coordinator_.barrier_.wait([this, &name] {
    auto now = std::chrono::steady_clock::now();
    coordinator_.durations_[std::move(name)] = now - coordinator_.start_;
  });
}

template <typename Iter, typename Work>
inline void Context::execute_synchronized_blockwise(Iter begin, Iter end,
                                                    Work work) {
  coordinator_.barrier_.wait([this] { coordinator_.reset_index(); });
  coordinator_.barrier_.wait();
  Iter block_begin = begin + coordinator_.get_next_index();
  while (block_begin < end) {
    Iter block_end = block_begin + coordinator_.block_size < end
                         ? block_begin + coordinator_.block_size
                         : end;
    work(block_begin, block_end);
    block_begin = begin + coordinator_.get_next_index();
  }
  coordinator_.barrier_.wait();
}

template <typename Iter, typename Work>
inline void Context::execute_synchronized_blockwise_timed(std::string name,
                                                          Iter begin, Iter end,
                                                          Work work) {
  coordinator_.barrier_.wait([this] { coordinator_.reset_index(); });
  coordinator_.barrier_.wait(
      [this] { coordinator_.start_ = std::chrono::steady_clock::now(); });
  Iter block_begin = begin + coordinator_.get_next_index();
  while (block_begin < end) {
    Iter block_end = block_begin + coordinator_.block_size < end
                         ? block_begin + coordinator_.block_size
                         : end;
    work(block_begin, block_end);
    block_begin = begin + coordinator_.get_next_index();
  }
  coordinator_.barrier_.wait([this, &name] {
    auto now = std::chrono::steady_clock::now();
    coordinator_.durations_[std::move(name)] = now - coordinator_.start_;
  });
}
}  // namespace thread_coordination

#endif
