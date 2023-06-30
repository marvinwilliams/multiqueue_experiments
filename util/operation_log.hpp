#pragma once

#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <ctime>
#include <filesystem>
#include <optional>
#include <ostream>
#include <vector>

namespace operation_log {

static_assert(sizeof(unsigned long) >= sizeof(std::uint64_t), "unsigned long is too small");
static constexpr std::uint8_t ElemIdBits = 56;

constexpr auto extract_thread_id(unsigned long packed_value) noexcept {
    return static_cast<int>(packed_value >> ElemIdBits);
}

constexpr auto extract_elem_id(unsigned long packed_value) noexcept {
    return static_cast<std::size_t>(packed_value & ((1UL << ElemIdBits) - 1));
}

constexpr unsigned long start_value(int thread_id) noexcept {
    assert(thread_id >= 0 && static_cast<unsigned long>(thread_id) < (1UL << ElemIdBits));
    return static_cast<unsigned long>(thread_id) << (ElemIdBits);
}

struct Push {
    std::uint64_t tick;
    unsigned long key;
    friend bool operator<(Push const& a, Push const& b) noexcept {
        return a.tick < b.tick;
    }
};

struct Pop {
    std::uint64_t tick;
    unsigned long data;
    friend bool operator<(Pop const& a, Pop const& b) noexcept {
        return a.tick < b.tick;
    }
};

struct OperationLog {
    std::vector<Push> pushes;
    std::vector<Pop> pops;
    long long failed_pops = 0;
};

template <typename PriorityQueue>
class LoggingHandle : public PriorityQueue::handle_type {
    using base_type = typename PriorityQueue::handle_type;
    unsigned long value_;
    OperationLog log_;

    std::uint64_t get_tick() noexcept {
        return static_cast<std::uint64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    }

   public:
    explicit LoggingHandle(int id, PriorityQueue& pq) : base_type(pq.get_handle()), value_(start_value(id)) {
    }

    void reserve_push_log(std::size_t size) {
        log_.pushes.resize(size);
        log_.pushes.clear();
    }

    void reserve_pop_log(std::size_t size) {
        log_.pops.resize(size);
        log_.pops.clear();
    }

    void push(typename PriorityQueue::key_type const& key) {
        base_type::push({key, value_});
        auto tick = get_tick();
        log_.pushes.push_back({tick, key});
        ++value_;
    }

    std::optional<typename PriorityQueue::value_type> try_pop() {
        auto tick = get_tick();
        auto retval = base_type::try_pop();
        if (!retval) {
            ++log_.failed_pops;
            return std::nullopt;
        }
        log_.pops.push_back({tick, retval->second});
        return *retval;
    }

    [[nodiscard]] OperationLog const& get_log() noexcept {
        return log_;
    }

    OperationLog extract_log() noexcept {
        OperationLog log{std::move(log_)};
        return log;
    }
};

struct Metrics {
    std::size_t rank_error;
    std::size_t delay;
};

bool verify_logs(std::vector<OperationLog> const& logs);
void write_logs(std::vector<OperationLog> const& log, std::ostream& out);
std::vector<Metrics> replay_logs(std::vector<OperationLog> const& logs);
void write_metrics(std::vector<Metrics> const& metrics, std::ostream& out);

}  // namespace operation_log
