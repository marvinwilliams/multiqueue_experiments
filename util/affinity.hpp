#pragma once

#include "threading.hpp"

#include <cstddef>

namespace affinity {

struct individual_cores {
    std::size_t stride = 1;
    std::size_t offset = 0;

    threading::thread_config operator()(int id) const {
        threading::thread_config cfg;
        cfg.cpu_set.set(offset + static_cast<std::size_t>(id) * stride);
        return cfg;
    }
};

struct same_core {
    std::size_t core = 1;

    threading::thread_config operator()(int /*unused*/) const {
        threading::thread_config cfg;
        cfg.cpu_set.set(core);
        return cfg;
    }
};

struct NUMA {
    int cores_per_node = 4;
    int num_nodes = 16;

    threading::thread_config operator()(int id) const {
        threading::thread_config cfg;
        int cpu = 0;
        if (id % (2 * cores_per_node) < cores_per_node) {
            cpu = (id / cores_per_node) * (cores_per_node / 2) + (id % cores_per_node);
        } else {
            id -= cores_per_node;
            cpu = num_nodes * cores_per_node + (id / cores_per_node) * cores_per_node / 2 + (id % cores_per_node);
        }
        cfg.cpu_set.set(static_cast<std::size_t>(cpu));
        return cfg;
    }
};

}  // namespace affinity
