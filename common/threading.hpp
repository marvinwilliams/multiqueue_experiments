#ifndef THREADING_HPP_YHEJJOUQ
#define THREADING_HPP_YHEJJOUQ

#include <sched.h>

#include <bitset>
#include <chrono>
#include <cstdio>
#include <functional>
#include <iostream>
#include <memory>
#include <optional>
#include <system_error>
#include <variant>

namespace threading {

namespace detail {

extern "C" void *trampoline(void *invoker);

struct invoker_base {
    invoker_base() = default;
    invoker_base(invoker_base const &) = delete;
    invoker_base(invoker_base &&) = delete;
    invoker_base &operator=(invoker_base const &) = delete;
    invoker_base &operator=(invoker_base &&) = delete;
    virtual ~invoker_base() = default;
    virtual void operator()() noexcept = 0;
};

template <typename Callable, typename... Args>
struct invoker : invoker_base {
    using ArgTuple = std::tuple<Args...>;
    Callable f;
    ArgTuple args;

    invoker(Callable func, std::tuple<Args...> arg_tuple) : f{std::move(func)}, args{std::move(arg_tuple)} {
    }

    template <size_t... I>
    void invoke(std::index_sequence<I...> /*unused*/) {
        std::invoke(std::move(f), std::get<I>(args)...);
    }

    void operator()() noexcept override {
        invoke(std::make_index_sequence<std::tuple_size_v<ArgTuple>>());
    }
};

}  // namespace detail

namespace scheduling {
struct Default {
    static constexpr int PolicyId = SCHED_OTHER;
    static constexpr int priority = 0;
};

struct Fifo {
    static constexpr int PolicyId = SCHED_FIFO;

    int priority;

    constexpr explicit Fifo(int prio = sched_get_priority_max(PolicyId)) noexcept : priority{prio} {
    }
};

struct RoundRobin {
    static constexpr int PolicyId = SCHED_RR;

    int priority;

    constexpr explicit RoundRobin(int prio = sched_get_priority_max(PolicyId)) noexcept : priority{prio} {
    }
};

using Policy = std::variant<scheduling::Default, scheduling::Fifo, scheduling::RoundRobin>;

}  // namespace scheduling

struct thread_config {
    static constexpr auto cpu_setsize = CPU_SETSIZE;
    std::bitset<cpu_setsize> cpu_set;
    std::optional<scheduling::Policy> scheduling;
    bool detached = false;
};

class pthread {
    std::optional<pthread_t> thread_handle_ = std::nullopt;
    bool detached_ = false;

    static void init_attr_scheduling(pthread_attr_t &attr, scheduling::Policy policy);
    static void init_attr(pthread_attr_t &attr, thread_config const &config);

    template <typename Callable, typename... Args>
    void start_thread(pthread_attr_t *attr, Callable &&f, Args &&...args) {
        static_assert(std::is_invocable_v<std::decay_t<Callable>, std::decay_t<Args>...>,
                      "pthread arguments must be invocable after conversion to rvalues");
        using invoker_type = detail::invoker<std::decay_t<Callable>, std::decay_t<Args>...>;
        auto invoker_ptr = std::unique_ptr<detail::invoker_base>{
            new invoker_type(std::forward<Callable>(f), {std::forward<Args>(args)...})};
        thread_handle_ = pthread_t{};
        if (int rc = pthread_create(&(*thread_handle_), attr, detail::trampoline, invoker_ptr.get())) {
            throw std::system_error{rc, std::system_category(), "Failed to create thread: "};
        }
        (void)invoker_ptr.release();
    }

   public:
    template <typename Callable, typename... Args>
    pthread(thread_config const &config, Callable &&f, Args &&...args) : detached_(config.detached) {
        pthread_attr_t attr;
        init_attr(attr, config);
        try {
            start_thread(&attr, std::forward<Callable>(f), std::forward<Args>(args)...);
            pthread_attr_destroy(&attr);
        } catch (std::system_error &e) {
            pthread_attr_destroy(&attr);
            throw;
        }
    }

    template <typename Callable, typename... Args,
              typename = std::enable_if_t<!std::is_same_v<std::decay_t<Callable>, thread_config>>>
    explicit pthread(Callable &&f, Args &&...args) {
        start_thread(nullptr, std::forward<Callable>(f), std::forward<Args>(args)...);
    }

    pthread() noexcept = default;

    pthread(pthread const &) = delete;

    pthread(pthread &&other) noexcept;

    pthread &operator=(pthread const &other) = delete;

    pthread &operator=(pthread &&other) noexcept;

    ~pthread() noexcept;

    [[nodiscard]] bool joinable() const;

    void detach();

    void set_policy(scheduling::Policy policy);

    void set_priority(int priority);

    void pin_to_core(size_t core);

    void set_affinity(std::bitset<CPU_SETSIZE> const &cpu_set);

    bool join();

    bool try_join();

    bool join_for(std::chrono::milliseconds ms);

    void cancel();
};

std::chrono::milliseconds get_thread_runtime();

}  // namespace threading

#endif
