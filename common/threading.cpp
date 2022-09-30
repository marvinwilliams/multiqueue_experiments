#include "threading.hpp"

#include <pthread.h>
#include <cassert>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>
#include <variant>

namespace threading {

extern "C" void *trampoline(void *invoker) {
    auto ptr = std::unique_ptr<detail::invoker_base>{static_cast<detail::invoker_base *>(invoker)};
    (*ptr)();
    return nullptr;
}

pthread::pthread(pthread &&other) noexcept {
    std::swap(thread_handle_, other.thread_handle_);
    std::swap(detached_, other.detached_);
}

pthread &pthread::operator=(pthread &&other) noexcept {
    if (thread_handle_) {
        std::terminate();
    }
    std::swap(thread_handle_, other.thread_handle_);
    std::swap(detached_, other.detached_);
    return *this;
}

pthread::~pthread() noexcept {
    if (joinable()) {
        std::terminate();
    }
}

void pthread::init_attr_scheduling(pthread_attr_t &attr, scheduling::Policy policy) {
    if (int rc = pthread_attr_setinheritsched(&attr, PTHREAD_EXPLICIT_SCHED)) {
        throw std::system_error{rc, std::system_category(), "Failed to set explicit scheduling for thread: "};
    }
    std::visit(
        [&attr](auto const &policy) {
            if (int rc = pthread_attr_setschedpolicy(&attr, policy.PolicyId)) {
                throw std::system_error{rc, std::system_category(), "Failed to set scheduling policy: "};
            }
            sched_param param{};
            param.sched_priority = policy.priority;
            if (int rc = pthread_attr_setschedparam(&attr, &param)) {
                throw std::system_error{rc, std::system_category(), "Failed to set scheduling parameter: "};
            }
        },
        policy);
}

void pthread::init_attr(pthread_attr_t &attr, thread_config const &config) {
    if (int rc = pthread_attr_init(&attr)) {
        throw std::system_error{rc, std::system_category(), "Failed to init pthread attribute: "};
    }
    if (config.detached) {
        if (int rc = pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED)) {
            throw std::system_error{rc, std::system_category(), "Failed to set to detached state: "};
        }
    }
    if (config.scheduling) {
        init_attr_scheduling(attr, *config.scheduling);
    }
    if (!config.cpu_set.none()) {
        cpu_set_t set{};
        CPU_ZERO(&set);
        for (int i = 0; i < config.cpu_set.size(); ++i) {
            if (config.cpu_set[i]) {
                CPU_SET(i, &set);
            }
        }
        if (int rc = pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set)) {
            throw std::system_error{rc, std::system_category(), "Failed to set thread affinity: "};
        }
    }
}

bool pthread::joinable() const {
    return thread_handle_.has_value() && !detached_;
}

void pthread::detach() {
    if (!thread_handle_) {
        throw std::runtime_error{"Thread object does not represent an active thread of execution"};
    }
    if (detached_) {
        throw std::runtime_error{"Thread is already detached"};
    }
    int rc = pthread_detach(*thread_handle_);
    if (rc != 0) {
        throw std::system_error{rc, std::system_category(), "Failed to detach thread: "};
    }
    detached_ = true;
}

void pthread::set_policy(scheduling::Policy policy) {
    if (!thread_handle_) {
        throw std::runtime_error{"Thread object does not represent an active thread of execution"};
    }
    std::visit(
        [this](auto const &p) {
            sched_param param{};
            param.sched_priority = p.priority;
            if (int rc = pthread_setschedparam(*thread_handle_, p.PolicyId, &param)) {
                throw std::system_error{rc, std::system_category(), "Failed to set scheduling parameter: "};
            }
        },
        policy);
}

void pthread::set_priority(int priority) {
    if (!thread_handle_) {
        throw std::runtime_error{"Thread object does not represent an active thread of execution"};
    }
    if (int rc = pthread_setschedprio(*thread_handle_, priority)) {
        throw std::system_error{rc, std::system_category(), "Failed to set thread priority: "};
    }
}

void pthread::pin_to_core(std::size_t core) {
    if (!thread_handle_) {
        throw std::runtime_error{"Thread object does not represent an active thread of execution"};
    }
    cpu_set_t set;
    CPU_ZERO(&set);
    CPU_SET(core, &set);
    if (int rc = pthread_setaffinity_np(*thread_handle_, sizeof(cpu_set_t), &set)) {
        throw std::system_error{rc, std::system_category(), "Pin to core failed: "};
    }
}

void pthread::set_affinity(std::bitset<CPU_SETSIZE> const &cpu_set) {
    if (!thread_handle_) {
        throw std::runtime_error{"Thread object does not represent an active thread of execution"};
    }
    cpu_set_t set;
    CPU_ZERO(&set);
    for (std::size_t i = 0; i < cpu_set.size(); ++i) {
        if (cpu_set[i]) {
            CPU_SET(i, &set);
        }
    }
    int rc = pthread_setaffinity_np(*thread_handle_, sizeof(cpu_set_t), &set);
    if (rc != 0) {
        throw std::system_error{rc, std::system_category(), "Failed to set thread affinity"};
    }
}

bool pthread::join() {
    if (!thread_handle_) {
        throw std::runtime_error{"Thread object does not represent an active thread of execution"};
    }
    if (detached_) {
        throw std::runtime_error{"Cannot join detached thread"};
    }
    void *retval = nullptr;
    if (int rc = pthread_join(*thread_handle_, &retval)) {
        throw std::system_error{rc, std::system_category(), "Failed to join thread: "};
    }
    thread_handle_ = std::nullopt;
    return retval != PTHREAD_CANCELED;
}

bool pthread::try_join() {
    if (!thread_handle_) {
        throw std::runtime_error{"Thread object does not represent an active thread of execution"};
    }
    if (detached_) {
        throw std::runtime_error{"Cannot join detached thread"};
    }
    void *retval = nullptr;
    int rc = pthread_tryjoin_np(*thread_handle_, &retval);
    if (rc == EBUSY) {
        return false;
    }
    if (rc != 0) {
        throw std::system_error{rc, std::system_category(), "Failed to try to join: "};
    }
    thread_handle_ = std::nullopt;
    return true;
}

bool pthread::join_for(std::chrono::milliseconds ms) {
    if (!thread_handle_) {
        throw std::runtime_error{"Thread object does not represent an active thread of execution"};
    }
    if (detached_) {
        throw std::runtime_error{"Cannot join detached thread"};
    }
    void *retval = nullptr;
    timespec ts{};
    if (clock_gettime(CLOCK_REALTIME, &ts) == -1) {
        throw std::system_error{errno, std::system_category(), "Failed to get time: "};
    }
    auto s = std::chrono::duration_cast<std::chrono::seconds>(ms);
    ts.tv_sec += s.count();
    ts.tv_nsec += std::chrono::nanoseconds{ms - s}.count();
    int rc = pthread_timedjoin_np(*thread_handle_, &retval, &ts);
    if (rc == ETIMEDOUT) {
        return false;
    }
    if (rc != 0) {
        throw std::system_error{rc, std::system_category(), "Failed to join timed: "};
    }
    thread_handle_ = std::nullopt;
    return true;
}

void pthread::cancel() {
    if (!thread_handle_) {
        throw std::runtime_error{"Thread object does not represent an active thread of execution"};
    }
    int rc = pthread_cancel(*thread_handle_);
    if (rc != 0) {
        throw std::system_error{rc, std::system_category(), "Failed to cancel thread: "};
    }
}

std::chrono::milliseconds get_thread_runtime() {
    clockid_t cid = 0;
    auto rc = pthread_getcpuclockid(pthread_self(), &cid);
    if (rc != 0) {
        throw std::system_error{rc, std::system_category(), "Failed to get thread clock id"};
    }
    timespec ts{};
    if (clock_gettime(cid, &ts) == -1) {
        throw std::system_error{errno, std::system_category(), "Failed to get thread runtime"};
    }
    return std::chrono::seconds{ts.tv_sec} +
        std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::nanoseconds{ts.tv_nsec});
}

}  // namespace threading
