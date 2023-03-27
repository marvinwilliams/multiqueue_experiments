#pragma once

#include <cassert>
#include <cstdint>
#include <optional>
#include <vector>

using key_type = unsigned long;
using value_type = unsigned long;
using tick_type = std::uint64_t;

namespace packed_value {

static constexpr unsigned int ThreadIdBits = 7;
static constexpr unsigned int OpIdBits = 32 - ThreadIdBits;  // Use maximum of 32 bits
static constexpr int MaxThreadId = 1 << ThreadIdBits;
static constexpr value_type MaxOpId = value_type{1} << OpIdBits;

constexpr value_type pack(int thread_id, std::size_t op_id) noexcept {
    assert(thread_id < MaxThreadId);
    assert(op_id < MaxOpId);
    return static_cast<value_type>(op_id) | (static_cast<value_type>(thread_id) << (OpIdBits));
}

constexpr int thread_id(value_type packed_value) noexcept {
    return static_cast<int>(packed_value >> OpIdBits);
}

constexpr std::size_t op_id(value_type packed_value) noexcept {
    return static_cast<std::size_t>(packed_value) & (MaxOpId - 1);
}

}  // namespace packed_value

struct PushLogEntry {
    tick_type tick;
    key_type key;
};

struct PopLogEntry {
    struct Payload {
        int from_thread_id;
        std::size_t op_id;
    };
    tick_type tick;
    Payload payload{};

    PopLogEntry(tick_type t, value_type v)
        : tick{t}, payload{Payload{packed_value::thread_id(v), packed_value::op_id(v)}} {
    }
};

using PushLogType = std::vector<std::vector<PushLogEntry>>;
using PopLogType = std::vector<std::vector<PopLogEntry>>;
