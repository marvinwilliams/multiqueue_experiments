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

struct Push {
    std::uint64_t tick;
    unsigned long key;
    friend bool operator<(Push const& a, Push const& b) noexcept {
        return a.tick < b.tick;
    }
};

struct Pop {
    std::uint64_t tick;
    std::size_t index;
    friend bool operator<(Pop const& a, Pop const& b) noexcept {
        return a.tick < b.tick;
    }
};

struct OperationLog {
    std::vector<Push> pushes;
    std::vector<Pop> pops;
};

inline std::uint64_t get_tick() noexcept {
    return static_cast<std::uint64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
}

struct Metrics {
    std::size_t rank_error;
    std::size_t delay;
};

bool verify_logs(OperationLog const& logs);
void write_logs(OperationLog const& log, std::ostream& out);
std::vector<Metrics> replay_logs(OperationLog logs);
void write_metrics(std::vector<Metrics> const& metrics, std::ostream& out);

}  // namespace operation_log
