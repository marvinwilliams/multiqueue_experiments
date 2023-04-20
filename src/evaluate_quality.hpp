#pragma once

#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <vector>

namespace quality {

using key_type = unsigned long;
using value_type = unsigned long;
using clock_type = std::chrono::steady_clock;
using tick_type = clock_type::time_point;

namespace packed_value {

static constexpr std::uint8_t ThreadIdBits = 7U;
static constexpr std::uint8_t ElemIdBits = 32U - ThreadIdBits;  // Use maximum of 32 bits
static constexpr std::uint32_t MaxThreadId = 1U << ThreadIdBits;
static constexpr std::uint32_t MaxElemId = 1U << ElemIdBits;

constexpr int get_thread_id(value_type packed_value) noexcept {
    return static_cast<int>(packed_value >> ElemIdBits);
}

constexpr std::size_t get_elem_id(value_type packed_value) noexcept {
    return static_cast<std::size_t>(packed_value) & (MaxElemId - 1);
}

constexpr value_type pack(int thread_id, std::size_t elem_id) noexcept {
    assert(thread_id >= 0 && thread_id < static_cast<int>(MaxThreadId));
    assert(elem_id < MaxElemId);
    std::uint32_t packed =
        static_cast<std::uint32_t>(elem_id) | (static_cast<std::uint32_t>(thread_id) << (ElemIdBits));
    assert(get_thread_id(packed) == thread_id && get_elem_id(packed) == elem_id);
    return static_cast<value_type>(packed);
}

}  // namespace packed_value

struct PushLogEntry {
    tick_type tick;
    key_type key;
};

struct PopLogEntry {
    tick_type tick;
    int thread_id;
    std::size_t elem_id;

    PopLogEntry(tick_type t, value_type v)
        : tick{t}, thread_id{packed_value::get_thread_id(v)}, elem_id{packed_value::get_elem_id(v)} {
    }
};

using PushLogType = std::vector<std::vector<PushLogEntry>>;
using PopLogType = std::vector<std::vector<PopLogEntry>>;

void write_logs(PushLogType const& push_log, PopLogType const& pop_log, std::filesystem::path const& file);
bool fix_and_verify_logs(PushLogType& push_log, PopLogType const& pop_log);
void write_histogram(PushLogType const& push_log, PopLogType const& pop_log, std::filesystem::path const& file);

}  // namespace quality
