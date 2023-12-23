#pragma once

#ifdef WITH_PAPI
#include <papi.h>
#endif

#include <ostream>

#ifdef CORES_PER_NUMA_NODE
static constexpr auto cores_per_numa_node = CORES_PER_NUMA_NODE;
#else
static constexpr auto cores_per_numa_node = 4;
#endif
#ifdef NUM_NUMA_NODES
static constexpr auto num_numa_nodes = NUM_NUMA_NODES;
#else
static constexpr auto num_numa_nodes = 16;
#endif

static std::ostream& write_build_info(std::ostream& out) {
    out << "Built on " << __DATE__ << ' ' << __TIME__ << " with\n";
#if defined(__clang__)
    out << "  Clang " << __clang_version__ << '\n';
#elif defined(__GNUC__)
    out << "  GCC " << __VERSION__ << '\n';
#else
    out << "  Unknown compiler\n";
#endif
#ifdef NDEBUG
    out << "  NDEBUG defined\n";
#else
    out << "  NDEBUG not defined\n";
#endif
#if defined WITH_PAPI
    out << "  PAPI " << PAPI_VER_CURRENT << '\n';
#else
    out << "  PAPI unsupported\n";
#endif
    return out;
}
