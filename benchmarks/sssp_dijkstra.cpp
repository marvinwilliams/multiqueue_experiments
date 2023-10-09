#include "build_info.hpp"
#include "graph.hpp"
#include "task.hpp"
#include "termination_detection.hpp"
#include "wrapper/priority.hpp"
#include "wrapper/selector.hpp"

#include "cxxopts.hpp"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <x86intrin.h>
#include <array>
#include <atomic>
#include <cassert>
#include <charconv>
#include <chrono>
#include <condition_variable>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <new>
#include <numeric>
#include <random>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

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
using pq_type = PriorityQueue<unsigned long, unsigned long, Priority::Min>;
using handle_type = pq_type::handle_type;

struct Settings {
    int num_threads = 4;
    std::filesystem::path graph_file;
    std::filesystem::path distance_file;
    int seed = 1;
    pq_type::config_type pq_settings{};
};

void write_settings(Settings const& settings, std::ostream& out) {
    out << "Threads: " << settings.num_threads << '\n'
        << "Graph: " << settings.graph_file.string() << '\n'
        << "Seed: " << settings.seed;
    out << "\n\n";
}

struct ThreadStats {
    std::pair<std::chrono::high_resolution_clock::time_point, std::chrono::high_resolution_clock::time_point> work_time;
    long long pushed_nodes{0};
    long long ignored_nodes{0};
    long long processed_nodes{0};
};

struct SharedData {
    struct alignas(L1_CACHE_LINESIZE) AtomicDistance {
        std::atomic<long long> value{std::numeric_limits<long long>::max()};
    };
    Graph graph;
    std::vector<AtomicDistance> shortest_distances;
    termination_detection::Data termination_detection_data{};

    SharedData(std::size_t num_nodes) : shortest_distances(num_nodes) {
    }

    bool update_distance(std::size_t index, long long current, long long target) noexcept {
        while (target < current) {
            if (shortest_distances[index].value.compare_exchange_weak(current, target, std::memory_order_relaxed)) {
                return true;
            }
        }
        return false;
    }
};

bool process_node(handle_type& handle, ThreadStats& stats, SharedData& data) {
    auto node = handle.try_pop();
    if (!node) {
        return false;
    }
    auto current_distance = data.shortest_distances[node->second].value.load(std::memory_order_relaxed);
    if (static_cast<long long>(node->first) > current_distance) {
        ++stats.ignored_nodes;
        return true;
    }
    ++stats.processed_nodes;
    for (auto i = data.graph.nodes[node->second]; i < data.graph.nodes[node->second + 1]; ++i) {
        auto target = data.graph.edges[i].target;
        auto d = static_cast<long long>(node->first) + data.graph.edges[i].weight;
        auto old_d = data.shortest_distances[target].value.load(std::memory_order_relaxed);
        if (data.update_distance(target, old_d, d)) {
            handle.push({d, target});
            ++stats.pushed_nodes;
        }
    }
    return true;
}

ThreadStats benchmark_thread(task::Control tc, pq_type& pq, SharedData& data) {
    ThreadStats stats;
    auto handle = pq.get_handle();
    if (tc.id() == 0) {
        data.shortest_distances[0].value = 0;
        handle.push({0, 0});
        ++stats.pushed_nodes;
    }
    tc.synchronize();
    auto start_time = std::chrono::high_resolution_clock::now();
    while (termination_detection::try_do(tc.num_threads(), data.termination_detection_data,
                                         [&]() { return process_node(handle, stats, data); })) {
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    tc.synchronize();
    stats.work_time = {start_time, end_time};
    return stats;
}

void write_stats(ThreadStats const& stats, std::ostream& out) {
    std::cout << "time,pushed,processed,ignored\n";
    out << (stats.work_time.second - stats.work_time.first).count() << ',' << stats.pushed_nodes << ','
        << stats.processed_nodes << ',' << stats.ignored_nodes << '\n';
}

bool verify_distances(SharedData const& shared_data) {
    for (std::size_t i = 0; i < shared_data.graph.num_nodes(); ++i) {
        for (std::size_t j = shared_data.graph.nodes[i]; j < shared_data.graph.nodes[i + 1]; ++j) {
            auto d = shared_data.shortest_distances[i].value + shared_data.graph.edges[j].weight;
            if (d < shared_data.shortest_distances[shared_data.graph.edges[j].target].value) {
                return false;
            }
        }
    }
    return true;
}

bool run_benchmark(Settings const& settings) {
    std::ofstream distance_out;
    if (!settings.distance_file.empty()) {
        distance_out = std::ofstream(settings.distance_file);
        if (!distance_out) {
            std::cerr << "Error: Could not open file " << settings.distance_file << " for writing" << std::endl;
            return false;
        }
    }

    std::clog << "Reading graph..." << std::endl;
    Graph graph;
    try {
        graph = Graph(settings.graph_file);
    } catch (std::runtime_error const& e) {
        std::clog << "Error: " << e.what() << std::endl;
        return false;
    }
    SharedData shared_data{graph.num_nodes()};
    shared_data.graph = std::move(graph);

    auto pq = pq_type(settings.num_threads, shared_data.graph.num_nodes(), settings.pq_settings);
    std::clog << "Priority queue: ";
    pq.describe(std::clog) << '\n';
    std::vector<ThreadStats> all_stats(static_cast<std::size_t>(settings.num_threads));
    affinity::NUMA numa_affinity{cores_per_numa_node, num_numa_nodes};
    std::clog << "Solving..." << std::endl;
    task::Runner runner{numa_affinity, settings.num_threads, [&](auto tc) {
                            all_stats[static_cast<std::size_t>(tc.id())] = benchmark_thread(tc, pq, shared_data);
                        }};
    runner.wait();

    if (distance_out.is_open()) {
        std::clog << "Writing distances..." << std::endl;
        for (std::size_t i = 0; i < shared_data.shortest_distances.size(); ++i) {
            distance_out << i << ' ' << shared_data.shortest_distances[i].value << '\n';
        }
        distance_out.close();
    }
    std::clog << "Finished\n" << std::endl;
    auto accum_stats =
        std::accumulate(all_stats.begin() + 1, all_stats.end(), all_stats.front(), [](auto accum, auto const& e) {
            accum.work_time.first = std::min(accum.work_time.first, e.work_time.first);
            accum.work_time.second = std::max(accum.work_time.second, e.work_time.second);
            accum.pushed_nodes += e.pushed_nodes;
            accum.processed_nodes += e.processed_nodes;
            accum.ignored_nodes += e.ignored_nodes;
            return accum;
        });
    std::clog << "Time (s): " << std::setprecision(3)
              << std::chrono::duration<double>(accum_stats.work_time.second - accum_stats.work_time.first).count()
              << '\n';
    std::clog << "Pushed nodes: " << accum_stats.pushed_nodes << '\n';
    std::clog << "Processed nodes: " << accum_stats.processed_nodes << '\n';
    std::clog << "Ignored nodes: " << accum_stats.ignored_nodes << '\n';
    if (accum_stats.processed_nodes + accum_stats.ignored_nodes != accum_stats.pushed_nodes) {
        std::clog << "Error: " << accum_stats.pushed_nodes - (accum_stats.processed_nodes + accum_stats.ignored_nodes)
                  << " node(s) were not popped" << std::endl;
        return false;
    }
    if (!verify_distances(shared_data)) {
        std::clog << "Error: Invalid distances" << std::endl;
        return false;
    }
    write_stats(accum_stats, std::cout);
    return true;
}

int main(int argc, char* argv[]) {
    write_build_info(std::clog);
    std::clog << "\nCommand line:";
    for (int i = 0; i < argc; ++i) {
        std::clog << ' ' << argv[i];
    }
    std::clog << '\n' << '\n';

    Settings settings{};
    cxxopts::Options cmd(argv[0]);
    // clang-format off
    cmd.add_options()
      ("j,threads", "The number of threads", cxxopts::value<int>(settings.num_threads), "NUMBER")
      ("file", "The input graph", cxxopts::value<std::filesystem::path>(settings.graph_file), "PATH")
      ("o,distance-file", "Path to write the distances to", cxxopts::value<std::filesystem::path>(settings.distance_file), "PATH")
      ("h,help", "Print this help");
    // clang-format on
    pq_type::add_options(cmd, settings.pq_settings);
    cmd.parse_positional({"file"});

    try {
        auto args = cmd.parse(argc, argv);
        if (args.count("help") > 0) {
            std::cerr << cmd.help() << std::endl;
            return EXIT_SUCCESS;
        }
    } catch (cxxopts::OptionParseException const& e) {
        std::cerr << "Error parsing arguments: " << e.what() << '\n';
        std::cerr << cmd.help() << std::endl;
        return EXIT_FAILURE;
    }

    write_settings(settings, std::clog);

    bool success = run_benchmark(settings);

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
