#pragma once

#include "thread_coordination.hpp"
#ifdef WITH_OPERATION_LOG
#include "operation_log.hpp"

#include <chrono>
#include <ctime>
#endif

#ifdef WITH_PAPI
#include <papi.h>
#endif

#include <sstream>
#include <string>

static std::string get_build_info() {
    std::stringstream ss;
    ss << "Built on " << __DATE__ << ' ' << __TIME__ << " with\n";
#if defined(__clang__)
    ss << "  Clang " << __clang_version__ << '\n';
#elif defined(__GNUC__)
    ss << "  GCC " << __VERSION__ << '\n';
#else
    ss << "  Unknown compiler\n";
#endif
#ifdef NDEBUG
    ss << "  NDEBUG defined\n";
#else
    ss << "  NDEBUG not defined\n";
#endif
#if defined WITH_PAPI
    ss << "  PAPI " << PAPI_VER_CURRENT << '\n';
#else
    ss << "  PAPI unsupported\n";
#endif
#ifdef WITH_OPERATION_LOG
    ss << " Operation logging enabled\n";
#endif
    return ss.str();
}

template <typename PriorityQueue>
class HandleWrapper {
   public:
    using pq_type = PriorityQueue;
    using key_type = typename pq_type::key_type;
    using value_type = typename pq_type::value_type;
    using handle_type = typename pq_type::handle_type;

   private:
    handle_type handle_;
#ifdef WITH_OPERATION_LOG
    operation_log::OperationLog log_;

    std::uint64_t get_tick() noexcept {
        timespec ts{};
        clock_gettime(CLOCK_MONOTONIC_RAW, &ts);
        return static_cast<std::uint64_t>(ts.tv_sec) * std::chrono::nanoseconds{std::chrono::seconds{1}}.count() +
            static_cast<std::uint64_t>(ts.tv_nsec);
    }
#endif

   public:
    BenchmarkContext(handle_type handle) : handle_(std::move(handle)) {
    }

    void push_pq(key_type const& key) {
#ifdef WITH_OPERATION_LOG
        auto value = operation_log::pack(ctx_.get_id(), log_.pushes.size());
        handle_.push({key, value});
        auto tick = get_tick();
        log_.pushes.push_back({tick, key});
#else
        handle_.push({key, key});
#endif
    }

    template <typename Iter>
    void push_pq(Iter begin, Iter end) {
        for (auto it = begin; it != end; ++it) {
            push_pq(*it);
        }
    }

    bool try_pop_pq(value_type& retval) {
#ifdef WITH_OPERATION_LOG
        auto tick = get_tick();
        if (!handle_.try_pop(retval)) {
            return false;
        }
        log_.pops.push_back({tick, retval.second});
        return true;
#else
        return handle_.try_pop(retval);
#endif
    }

    bool try_pop_pq() {
        value_type retval;
        return try_pop_pq(retval);
    }

    value_type pop_pq() {
        value_type retval;
        while (!try_pop_pq(retval)) {
        }
        return retval;
    }

    auto get_pq_stats() const noexcept {
        return handle_.get_stats();
    }

#ifdef WITH_OPERATION_LOG
    operation_log::OperationLog& log() noexcept {
        return log_;
    }
#endif
};
