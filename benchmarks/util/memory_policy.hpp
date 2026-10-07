#pragma once

#include "util/thread_coordination.hpp"

#include <linux/mempolicy.h>
#include <sys/syscall.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <climits>
#include <cstdlib>
#include <cstddef>
#include <fstream>
#include <set>
#include <string>
#include <vector>

namespace memory_policy {

enum Policy : int {
    FirstTouch = 0,
    InterleaveUsed = 1,
    InterleaveAll = 2,
};

inline constexpr int max_id = InterleaveAll;
inline constexpr int default_id = InterleaveUsed;

inline char const* name(int id) noexcept {
    switch (id) {
        case FirstTouch:
            return "First touch";
        case InterleaveUsed:
            return "Interleave over used nodes";
        case InterleaveAll:
            return "Interleave over all nodes";
        default:
            return "Unknown";
    }
}

namespace detail {

inline std::vector<std::size_t> parse_list(std::string const& list) {
    std::vector<std::size_t> ids;
    std::size_t pos = 0;
    while (pos < list.size()) {
        auto end = list.find(',', pos);
        if (end == std::string::npos) {
            end = list.size();
        }
        auto range = list.substr(pos, end - pos);
        if (!range.empty()) {
            auto dash = range.find('-');
            auto first = std::strtoul(range.c_str(), nullptr, 10);
            auto last = dash == std::string::npos ? first : std::strtoul(range.c_str() + dash + 1, nullptr, 10);
            for (auto i = first; i <= last; ++i) {
                ids.push_back(i);
            }
        }
        pos = end + 1;
    }
    return ids;
}

inline std::vector<std::size_t> nodes_with_memory() {
    std::string list;
    std::ifstream("/sys/devices/system/node/has_memory") >> list;
    auto nodes = parse_list(list);
    if (nodes.empty()) {
        nodes.push_back(0);
    }
    return nodes;
}

}  // namespace detail

inline std::vector<std::size_t> nodes(int policy, int affinity, int num_threads) {
    if (policy == FirstTouch) {
        return {};
    }
    auto with_memory = detail::nodes_with_memory();
    if (policy == InterleaveAll) {
        return with_memory;
    }
    auto cpus = thread_coordination::affinity::cpu_assignment(affinity, num_threads);
    if (cpus.empty()) {
        return with_memory;
    }
    std::set<std::size_t> used;
    for (auto cpu : cpus) {
        used.insert(thread_coordination::affinity::numa_node(cpu));
    }
    std::vector<std::size_t> result;
    std::set_intersection(used.begin(), used.end(), with_memory.begin(), with_memory.end(),
                          std::back_inserter(result));
    if (result.empty()) {
        return with_memory;
    }
    return result;
}

[[nodiscard]] inline int apply(std::vector<std::size_t> const& nodes) {
    if (nodes.empty()) {
        return 0;
    }
    constexpr std::size_t bits = sizeof(unsigned long) * CHAR_BIT;
    std::array<unsigned long, 16> mask{};
    for (auto node : nodes) {
        if (node >= mask.size() * bits) {
            return EINVAL;
        }
        mask[node / bits] |= 1UL << (node % bits);
    }
    if (syscall(SYS_set_mempolicy, MPOL_INTERLEAVE, mask.data(), mask.size() * bits) != 0) {
        return errno;
    }
    return 0;
}

}  // namespace memory_policy
