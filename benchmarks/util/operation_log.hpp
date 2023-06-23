#pragma once

#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <ostream>
#include <vector>

namespace operation_log {

static_assert(sizeof(unsigned long) >= sizeof(std::uint64_t), "unsigned long is too small");
static constexpr std::uint8_t ElemIdBits = 56;

constexpr int extract_thread_id(unsigned long packed_value) noexcept {
    return static_cast<int>(packed_value >> ElemIdBits);
}

constexpr std::size_t extract_elem_id(unsigned long packed_value) noexcept {
    return static_cast<std::size_t>(packed_value & ((1UL << ElemIdBits) - 1));
}

constexpr unsigned long pack(int thread_id, std::size_t elem_id) noexcept {
    assert(thread_id >= 0 && static_cast<unsigned long>(thread_id) < (1UL << ElemIdBits));
    assert(elem_id < (1UL << ElemIdBits));
    unsigned long packed = elem_id | (static_cast<unsigned long>(thread_id) << (ElemIdBits));
    assert(extract_thread_id(packed) == thread_id && extract_elem_id(packed) == elem_id);
    return packed;
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
};

struct Histogram {
    std::vector<std::size_t> ranks;
    std::vector<std::size_t> delays;
};

void verify_logs(std::vector<OperationLog> const& logs);
void write_logs(std::vector<OperationLog> const& log, std::ostream& out);
Histogram to_histogram(std::vector<OperationLog> const& logs);
void write_histogram(Histogram const& histogram, std::ostream& out);

}  // namespace operation_log
