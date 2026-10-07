#pragma once

#include "util/json.hpp"

#ifdef WITH_JEMALLOC
#include <jemalloc/jemalloc.h>
#endif
#include <sys/resource.h>
#include <unistd.h>

#ifdef WITH_JEMALLOC
#include <cstddef>
#include <cstdint>
#endif
#include <fstream>
#include <initializer_list>
#include <utility>

namespace memory_stats {

inline long long current_rss_bytes() {
    long long pages = 0;
    long long resident = -1;
    std::ifstream("/proc/self/statm") >> pages >> resident;
    return resident < 0 ? -1 : resident * sysconf(_SC_PAGESIZE);
}

inline long long peak_rss_bytes() {
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0) {
        return -1;
    }
    return static_cast<long long>(usage.ru_maxrss) * 1024;
}

inline long long allocated_bytes() {
#ifdef WITH_JEMALLOC
    std::uint64_t epoch = 1;
    std::size_t size = sizeof(epoch);
    if (mallctl("epoch", &epoch, &size, &epoch, size) != 0) {
        return -1;
    }
    std::size_t allocated = 0;
    size = sizeof(allocated);
    if (mallctl("stats.allocated", &allocated, &size, nullptr, 0) != 0) {
        return -1;
    }
    return static_cast<long long>(allocated);
#else
    return -1;
#endif
}

struct Snapshot {
    long long rss_bytes = -1;
    long long allocated_bytes = -1;

    static Snapshot take() {
        return {current_rss_bytes(), memory_stats::allocated_bytes()};
    }

    void write_json(json::Object& obj) const {
        obj.entry("rss_bytes", rss_bytes);
        if (allocated_bytes >= 0) {
            obj.entry("allocated_bytes", allocated_bytes);
        }
    }
};

inline void write_json(json::Object& obj, std::initializer_list<std::pair<char const*, Snapshot>> phases) {
    for (auto const& [name, snapshot] : phases) {
        obj.object(name, [&s = snapshot](json::Object& o) { s.write_json(o); });
    }
    obj.entry("peak_rss_bytes", peak_rss_bytes());
}

}  // namespace memory_stats
