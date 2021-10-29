#ifndef THREAD_COORDINATION_HPP_7F7173D6
#define THREAD_COORDINATION_HPP_7F7173D6

#include "barrier.hpp"
#include "threading.hpp"

#include <pthread.h>

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace thread_coordination {

class ThreadCoordinator;

class alignas(L1_CACHE_LINESIZE) Context {
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

  template <typename StartFunc, typename EndFunc, typename Work>
  void execute_synchronized(StartFunc&& start_func, EndFunc&& end_func,
                            Work work);

  template <typename Iter, typename Work>
  void execute_synchronized_blockwise(Iter begin, Iter end,
                                      std::atomic<Iter>& counter, Work work);

  template <typename StartFunc, typename EndFunc, typename Iter, typename Work>
  void execute_synchronized_blockwise(StartFunc&& start_func,
                                      EndFunc&& end_func, Iter begin, Iter end,
                                      std::atomic<Iter>& counter, Work work);
};

class ThreadCoordinator {
  friend Context;
  unsigned int num_threads_;
  utils::barrier barrier_;
  std::vector<threading::pthread> threads_;

 public:
  explicit ThreadCoordinator(unsigned int num_threads)
      : num_threads_{num_threads}, barrier_{num_threads} {}

  unsigned int get_num_threads() const noexcept { return num_threads_; }

  template <typename Task, typename... Args>
  void run_task(Args&&... args) {
    threads_.clear();
    threads_.reserve(num_threads_);
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
};

inline unsigned int Context::get_num_threads() const noexcept {
  return coordinator_.get_num_threads();
}

inline void Context::synchronize() { coordinator_.barrier_.wait(); }

template <typename CompletionFunc>
inline void Context::synchronize(CompletionFunc&& f) {
  coordinator_.barrier_.wait(std::forward<CompletionFunc>(f));
}

template <typename StartFunc, typename EndFunc, typename Work>
inline void Context::execute_synchronized(StartFunc&& start_func,
                                          EndFunc&& end_func, Work work) {
  coordinator_.barrier_.wait(std::forward<StartFunc>(start_func));
  work();
  coordinator_.barrier_.wait(std::forward<EndFunc>(end_func));
}

template <typename Iter, typename Work>
inline void Context::execute_synchronized_blockwise(Iter begin, Iter end,
                                                    std::atomic<Iter>& counter,
                                                    Work work) {
  static constexpr std::ptrdiff_t block_size = 4096;
  coordinator_.barrier_.wait([begin, &counter] { counter = begin; });
  coordinator_.barrier_.wait();
  Iter current_begin, current_end;
  current_begin = counter.fetch_add(block_size, std::memory_order_relaxed);
  while (current_begin < end) {
    current_end =
        current_begin + std::min(block_size, std::distance(current_begin, end));
    work(current_begin, current_end);
    current_begin = counter.fetch_add(block_size, std::memory_order_relaxed);
  }
  coordinator_.barrier_.wait();
}

template <typename StartFunc, typename EndFunc, typename Iter, typename Work>
inline void Context::execute_synchronized_blockwise(StartFunc&& start_func,
                                                    EndFunc&& end_func,
                                                    Iter begin, Iter end,
                                                    std::atomic<Iter>& counter,
                                                    Work work) {
  static constexpr std::ptrdiff_t block_size = 4096;
  coordinator_.barrier_.wait([begin, &counter] { counter = begin; });
  coordinator_.barrier_.wait(std::forward<StartFunc>(start_func));
  Iter current_begin, current_end;
  current_begin = counter.fetch_add(block_size, std::memory_order_relaxed);
  while (current_begin < end) {
    current_end =
        current_begin + std::min(block_size, std::distance(current_begin, end));
    work(current_begin, current_end);
    current_begin = counter.fetch_add(block_size, std::memory_order_relaxed);
  }
  coordinator_.barrier_.wait(std::forward<EndFunc>(end_func));
}

}  // namespace thread_coordination

#endif
