#ifndef THREAD_COORDINATION_HPP_7F7173D6
#define THREAD_COORDINATION_HPP_7F7173D6

#include "threading.hpp"

#include <pthread.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <condition_variable>
#include <mutex>
#include <optional>
#include <system_error>
#include <utility>
#include <vector>

namespace thread_coordination {

class ThreadCoordinator;

class alignas(L1_CACHE_LINESIZE) Context {
  friend ThreadCoordinator;
  ThreadCoordinator& coordinator_;
  unsigned int id_;
  unsigned int num_threads_;
  // The next step to wait on
  unsigned int step_ = 0;

 private:
  Context(ThreadCoordinator& coordinator, unsigned int id,
          unsigned int num_threads)
      : coordinator_{coordinator}, id_{id}, num_threads_{num_threads} {}

 public:
  bool is_main() const noexcept { return id_ == 0; }

  unsigned int get_id() const noexcept { return id_; }

  unsigned int get_num_threads() const noexcept { return num_threads_; }

  void synchronize(unsigned int step);

  template <typename Callable>
  void synchronize(unsigned int step, Callable f);

  template <typename Callable, typename... Args>
  void execute_synchronized(unsigned int step, Callable f, Args&&... args);
};

class ThreadCoordinator {
  friend Context;
  unsigned int num_threads_;
  threading::barrier barrier_;
  std::vector<threading::pthread> threads_;

 public:
  explicit ThreadCoordinator(unsigned int num_threads)
      : num_threads_{num_threads}, barrier_{num_threads} {}

  template <typename Task, typename... Args>
  void run_task(Args&&... args) {
    threads_.clear();
    threads_.reserve(num_threads_);
    for (unsigned int i = 0; i < num_threads_; ++i) {
      Context ctx{*this, i, num_threads_};
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

inline void Context::synchronize(unsigned int step) {
  assert(step_ <= step);
  do {
    // One synchronization step consists of two barriers
    coordinator_.barrier_.wait();
    coordinator_.barrier_.wait();
  } while (step_++ != step);
}

template <typename CompletionFunc>
inline void Context::synchronize(unsigned int step, Callable&& f) {
  assert(step_ <= step);
  while (step_++ != step) {
    coordinator_.barrier_.wait();
  }
  coordinator_.barrier_.wait(std::forward<Callable>(f));
}

template <typename SetupFunc, typename TeardownFunc, typename Task,
          typename... Args>
inline void Context::execute_synchronized(unsigned int step, SetupFunc setup,
                                          TeardownFunc teardown, Task f,
                                          Args&&... args) {
  assert(step_ != step);
  if (step_ != step) {
    synchronize(step - 1);
  }
  ++step_;
  if (is_main()) {
    setup();
  }
  coordinator_.barrier_.wait();
  f(std::forward<Args>(args)...);
  coordinator_.barrier_.wait();
  if (is_main()) {
    setup();
  }
}

}  // namespace thread_coordination

#endif
