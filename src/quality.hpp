#pragma once

#include "priority_queue_factory.hpp"

#include <cassert>
#include <cstdint>

using key_type = unsigned long;
using value_type = unsigned long;
using tick_type = std::uint64_t;

using PriorityQueue = util::priority_queue_type<key_type, value_type, true>;
using Handle = typename PriorityQueue::handle_type;

namespace packed_value {

static constexpr unsigned int ThreadIdBits = 7;
static constexpr unsigned int OpIdBits = 32 - ThreadIdBits; // Use maximum of 32 bits
static constexpr unsigned int MaxThreadId = 1UL << ThreadIdBits;
static constexpr value_type MaxOpId = value_type{1} << OpIdBits;

static constexpr value_type pack(unsigned int thread_id, std::size_t op_id) noexcept {
    assert(thread_id < MaxThreadId);
    assert(op_id < MaxOpId);
    return value_type{op_id} | (value_type{thread_id} << (OpIdBits));
}

static constexpr unsigned int thread_id(value_type packed_value) noexcept {
    return static_cast<unsigned int>(packed_value >> OpIdBits);
}

static constexpr std::size_t op_id(value_type packed_value) noexcept {
    return static_cast<std::size_t>(packed_value) & (MaxOpId - 1);
}

}  // namespace packed_value

struct PushLogEntry {
    tick_type tick;
    key_type key;
};

struct PopLogEntry {
    struct Payload {
        unsigned int from_thread_id;
        std::size_t op_id;
    };
    tick_type tick;
    std::optional<Payload> payload{};

    PopLogEntry(tick_type t, value_type v)
        : tick{t}, payload{Payload{packed_value::thread_id(v), packed_value::op_id(v)}} {
    }

    PopLogEntry(tick_type t) : tick{t} {
    }
};

using PushLogType = std::vector<std::vector<PushLogEntry>>;
using PopLogType = std::vector<std::vector<PopLogEntry>>;
