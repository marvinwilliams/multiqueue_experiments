#pragma once

#include "util/build_info.hpp"
#include "util/json.hpp"
#include "util/thread_coordination.hpp"

#include <cxxopts.hpp>

#include <iostream>
#include <ostream>

namespace benchmark {

struct BaseSettings {
    int num_threads = 4;
    int affinity = 6;

    void write_human_readable(std::ostream& out) const {
        out << "Threads: " << num_threads << '\n';
        out << "Affinity: " << thread_coordination::affinity::name(affinity) << '\n';
    }

    void write_json(json::Object& obj) const {
        obj.entry("num_threads", num_threads);
        obj.entry("affinity", affinity);
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
                "6: Far L1 Close L3)"
                , cxxopts::value<int>(affinity), "NUMBER");
        // clang-format on
    }
};

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

}  // namespace benchmark
