#pragma once

#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <vector>

namespace quality {

namespace packed_value {

static_assert(sizeof(unsigned long) >= sizeof(std::uint64_t), "unsigned long is too small");
static constexpr std::uint8_t ElemIdBits = 56;

constexpr int get_thread_id(unsigned long packed_value) noexcept {
    return static_cast<int>(packed_value >> ElemIdBits);
}

constexpr std::size_t get_elem_id(unsigned long packed_value) noexcept {
    return static_cast<std::size_t>(packed_value & ((1UL << ElemIdBits) - 1));
}

constexpr unsigned long pack(int thread_id, std::size_t elem_id) noexcept {
    assert(thread_id >= 0 && static_cast<unsigned long>(thread_id) < (1UL << ElemIdBits));
    assert(elem_id < (1UL << ElemIdBits));
    unsigned long packed = elem_id | (static_cast<unsigned long>(thread_id) << (ElemIdBits));
    assert(get_thread_id(packed) == thread_id && get_elem_id(packed) == elem_id);
    return packed;
}

}  // namespace packed_value

struct OperationLog {
    struct OperationLogEntry {
        std::uint64_t tick;
        unsigned long payload;
    };

    std::vector<OperationLogEntry> push_log;
    std::vector<OperationLogEntry> pop_log;
};

void write_logs(std::vector<OperationLog> const& logs, std::filesystem::path const& file);
void check_logs(std::vector<OperationLog> const& logs);
void write_histogram(std::vector<OperationLog> const& log, std::filesystem::path const& file);

}  // namespace quality
