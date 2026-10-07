#pragma once

#include "util/build_info.hpp"
#include "util/json.hpp"
#include "util/memory_policy.hpp"
#include "util/memory_stats.hpp"
#include "util/thread_coordination.hpp"

#include <cxxopts.hpp>

#ifdef WITH_PAPI
#include <papi.h>
#include <pthread.h>
#endif

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace benchmark {

using clock_type = std::chrono::steady_clock;

struct Interval {
    clock_type::time_point start;
    clock_type::time_point end;
};

inline Interval span(std::vector<Interval> const& intervals) {
    Interval span = intervals.front();
    for (auto const& i : intervals) {
        span.start = std::min(span.start, i.start);
        span.end = std::max(span.end, i.end);
    }
    return span;
}

inline double seconds(std::vector<Interval> const& intervals) {
    auto s = span(intervals);
    return std::chrono::duration<double>(s.end - s.start).count();
}

inline long long nanoseconds(clock_type::duration d) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(d).count();
}

inline long long time_ns(std::vector<Interval> const& intervals) {
    auto s = span(intervals);
    return nanoseconds(s.end - s.start);
}

inline void write_interval(json::Object& obj, std::string const& prefix, Interval const& interval,
                           clock_type::time_point origin) {
    obj.entry(prefix + "start_time", nanoseconds(interval.start - origin));
    obj.entry(prefix + "end_time", nanoseconds(interval.end - origin));
}

struct BaseSettings {
    int num_threads = 4;
    int affinity = thread_coordination::affinity::default_id;
    int memory = memory_policy::default_id;

    [[nodiscard]] std::vector<std::size_t> memory_nodes() const {
        return memory_policy::nodes(memory, affinity, num_threads);
    }

    void apply_memory_policy() const {
        auto nodes = memory_nodes();
        if (int err = memory_policy::apply(nodes); err != 0) {
            throw std::system_error{err, std::system_category(), "Failed to set memory policy"};
        }
        if (!nodes.empty()) {
            std::clog << "Interleaving memory over NUMA nodes";
            for (auto node : nodes) {
                std::clog << ' ' << node;
            }
            std::clog << '\n';
        }
    }

    void write_human_readable(std::ostream& out) const {
        out << "Threads: " << num_threads << '\n';
        out << "Affinity: " << thread_coordination::affinity::name(affinity) << '\n';
        out << "Memory: " << memory_policy::name(memory) << '\n';
    }

    void write_json(json::Object& obj) const {
        obj.entry("num_threads", num_threads);
        obj.entry("affinity", affinity);
        obj.entry("memory_policy", memory);
        obj.array("memory_nodes", memory_nodes());
    }

    [[nodiscard]] bool validate() const {
        if (num_threads <= 0) {
            std::cerr << "Error: Number of threads must be greater than 0\n";
            return false;
        }
        if (affinity < 0 || affinity > thread_coordination::affinity::max_id) {
            std::cerr << "Error: Invalid affinity\n";
            return false;
        }
        if (memory < 0 || memory > memory_policy::max_id) {
            std::cerr << "Error: Invalid memory policy\n";
            return false;
        }
        if (thread_coordination::affinity::needs_distinct_cpus(affinity)) {
            auto num_cpus = thread_coordination::affinity::available_cpus().size();
            if (static_cast<std::size_t>(num_threads) > num_cpus) {
                std::cerr << "Error: Affinity '" << thread_coordination::affinity::name(affinity)
                          << "' pins each thread to its own CPU, but only " << num_cpus << " CPUs are available\n";
                return false;
            }
        }
        return true;
    }

    void register_cmd_options(cxxopts::Options& cmd) {
        // clang-format off
        cmd.add_options()
            ("j,threads", "Number of threads", cxxopts::value<int>(num_threads), "NUMBER")
            ("a,affinity", "CPU affinity ("
                "0: None, "
                "1: Thread Id, "
                "2: Same, "
                "3: Close caches, "
                "4: Far caches, "
                "5: Close L3 Far L1, "
                "6: Far L1 Close L3, "
                "7: Cores first (L3, NUMA, SMT))"
                , cxxopts::value<int>(affinity), "NUMBER")
            ("memory", "Memory placement ("
                "0: First touch, "
                "1: Interleave over the NUMA nodes of the pinned CPUs, "
                "2: Interleave over all NUMA nodes)"
                , cxxopts::value<int>(memory), "NUMBER");
        // clang-format on
    }
};

#ifdef WITH_PAPI
class Papi {
    [[noreturn]] static void abort_with(std::string const& message) {
        std::cerr << "Error: " << message << std::endl;
        std::abort();
    }

   public:
    std::vector<std::string> events;

    void register_cmd_options(cxxopts::Options& cmd) {
        cmd.add_options()("r,count-event", "Papi event to count", cxxopts::value<std::vector<std::string>>(events),
                          "STRING");
    }

    [[nodiscard]] bool validate() const {
        if (events.empty()) {
            return true;
        }
        if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) {
            std::cerr << "Error: Failed to initialize PAPI library\n";
            return false;
        }
        if (PAPI_thread_init(pthread_self) != PAPI_OK) {
            std::cerr << "Error: Failed to initialize PAPI thread support\n";
            return false;
        }
        for (auto const& name : events) {
            if (PAPI_query_named_event(name.c_str()) != PAPI_OK) {
                std::cerr << "Error: PAPI event '" << name << "' not available\n";
                return false;
            }
        }
        int event_set = PAPI_NULL;
        if (PAPI_create_eventset(&event_set) != PAPI_OK) {
            std::cerr << "Error: Failed to create PAPI event set\n";
            return false;
        }
        for (auto const& name : events) {
            int event = PAPI_NULL;
            if (PAPI_event_name_to_code(name.c_str(), &event) != PAPI_OK ||
                PAPI_add_event(event_set, event) != PAPI_OK) {
                std::cerr << "Error: PAPI events cannot be counted together ('" << name << "' failed)\n";
                return false;
            }
        }
        PAPI_cleanup_eventset(event_set);
        PAPI_destroy_eventset(&event_set);
        return true;
    }

    void write_human_readable(std::ostream& out) const {
        out << "PAPI events:";
        if (events.empty()) {
            out << " None";
        }
        for (auto const& e : events) {
            out << ' ' << e;
        }
        out << '\n';
    }

    void write_json(json::Object& obj) const {
        obj.array("papi_events", events);
    }

    [[nodiscard]] int create_event_set() const {
        int event_set = -1;
        if (events.empty()) {
            return event_set;
        }
        if (PAPI_register_thread() != PAPI_OK) {
            abort_with("Failed to register thread for PAPI");
        }
        event_set = PAPI_NULL;
        if (PAPI_create_eventset(&event_set) != PAPI_OK) {
            abort_with("Failed to create PAPI event set");
        }
        for (auto const& name : events) {
            int event = PAPI_NULL;
            if (PAPI_event_name_to_code(name.c_str(), &event) != PAPI_OK ||
                PAPI_add_event(event_set, event) != PAPI_OK) {
                abort_with("Failed to add PAPI event '" + name + '\'');
            }
        }
        return event_set;
    }

    void start(int event_set) const {
        if (!events.empty() && PAPI_start(event_set) != PAPI_OK) {
            abort_with("Failed to start performance counters");
        }
    }

    void stop(int event_set, std::vector<long long>& counters) const {
        if (events.empty()) {
            return;
        }
        counters.resize(events.size());
        if (PAPI_stop(event_set, counters.data()) != PAPI_OK) {
            abort_with("Failed to stop performance counters");
        }
    }

    void accumulate(std::vector<long long>& total, std::vector<long long> const& counters) const {
        total.resize(events.size());
        for (std::size_t i = 0; i < counters.size() && i < total.size(); ++i) {
            total[i] += counters[i];
        }
    }

    // Writes the counters as an object {event: count}; nothing if no events are counted
    void write_counters(json::Object& obj, std::string const& name, std::vector<long long> const& counters) const {
        if (events.empty()) {
            return;
        }
        obj.object(name, [&](json::Object& o) {
            for (std::size_t i = 0; i < events.size() && i < counters.size(); ++i) {
                o.entry(events[i], counters[i]);
            }
        });
    }
};
#else
class Papi {
   public:
    void register_cmd_options(cxxopts::Options& /*cmd*/) {
    }

    [[nodiscard]] bool validate() const {
        return true;
    }

    void write_human_readable(std::ostream& /*out*/) const {
    }

    void write_json(json::Object& /*obj*/) const {
    }

    [[nodiscard]] int create_event_set() const {
        return -1;
    }

    void start(int /*event_set*/) const {
    }

    void stop(int /*event_set*/, std::vector<long long>& /*counters*/) const {
    }

    void accumulate(std::vector<long long>& /*total*/, std::vector<long long> const& /*counters*/) const {
    }

    void write_counters(json::Object& /*obj*/, std::string const& /*name*/,
                        std::vector<long long> const& /*counters*/) const {
    }
};

#endif

template <typename PQ>
void write_header(int argc, char* argv[], std::ostream& out) {
    write_build_info(out);
    out << '\n';
    out << "= Priority queue =\n";
    PQ::write_human_readable(out);
    out << '\n';
    out << "= Command line =\n";
    for (int i = 0; i < argc; ++i) {
        out << argv[i];
        if (i != argc - 1) {
            out << ' ';
        }
    }
    out << '\n' << '\n';
}

template <typename PQ, typename Settings, typename Benchmark>
int run(int argc, char* argv[], Benchmark&& b) {
    write_header<PQ>(argc, argv, std::clog);

    cxxopts::Options cmd(argv[0]);
    cmd.add_options()("h,help", "Print this help");
    Settings settings{};
    settings.register_cmd_options(cmd);
    try {
        if (cmd.parse(argc, argv).count("help") > 0) {
            std::clog << cmd.help() << '\n';
            return EXIT_SUCCESS;
        }
    } catch (cxxopts::OptionException const& e) {
        std::cerr << "Error parsing command line: " << e.what() << '\n';
        std::cerr << "Use --help for usage information\n";
        return EXIT_FAILURE;
    }

    std::clog << "= Settings =\n";
    settings.write_human_readable(std::clog);
    std::clog << '\n';
    if (!settings.validate()) {
        return EXIT_FAILURE;
    }
    try {
        settings.base_settings.apply_memory_policy();
    } catch (std::system_error const& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
    std::clog << "= Running benchmark =\n";
    std::forward<Benchmark>(b)(settings);
    return EXIT_SUCCESS;
}

}  // namespace benchmark
