#pragma once

#include <cassert>
#include <chrono>
#include <cstdint>
#include <optional>
#include <vector>

using key_type = unsigned long;
using value_type = unsigned long;
using tick_type = std::chrono::steady_clock::time_point;

namespace packed_value {

static constexpr unsigned int ThreadIdBits = 7;
static constexpr unsigned int ElemIdBits = 32 - ThreadIdBits;  // Use maximum of 32 bits
static constexpr int MaxThreadId = 1 << ThreadIdBits;
static constexpr value_type MaxElemId = value_type{1} << ElemIdBits;

constexpr value_type pack(int thread_id, std::size_t elem_id) noexcept {
    assert(thread_id < MaxThreadId);
    assert(elem_id < MaxElemId);
    return static_cast<value_type>(elem_id) | (static_cast<value_type>(thread_id) << (ElemIdBits));
}

constexpr int thread_id(value_type packed_value) noexcept {
    return static_cast<int>(packed_value >> ElemIdBits);
}

constexpr std::size_t elem_id(value_type packed_value) noexcept {
    return static_cast<std::size_t>(packed_value) & (MaxElemId - 1);
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
        : tick{t}, thread_id{packed_value::thread_id(v)}, elem_id{packed_value::elem_id(v)} {
    }
};

using PushLogType = std::vector<std::vector<PushLogEntry>>;
using PopLogType = std::vector<std::vector<PopLogEntry>>;
