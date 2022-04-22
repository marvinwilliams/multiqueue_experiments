#include "utils/priority_queue_factory.hpp"
#include "utils/thread_coordination.hpp"
#include "utils/threading.hpp"
#include "utils/xoroshiro256starstar.hpp"

#include "cxxopts.hpp"

#include <x86intrin.h>
#include <array>
#include <atomic>
#include <cassert>
#include <chrono>
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

using distance_type = unsigned long;
using node_type = unsigned long;

struct Graph {
    struct Edge {
        node_type target;
        distance_type weight;
    };
    std::vector<node_type> nodes;
    std::vector<Edge> edges;
};

struct alignas(2 * L1_CACHE_LINESIZE) Distance {
    std::atomic_ulong value;
};

using PriorityQueue =
    typename util::PriorityQueueFactory<distance_type, node_type>::type;

struct Settings {
    std::filesystem::path graph_file;
    std::filesystem::path output;
    unsigned int num_threads = 4;
    std::uint64_t seed = 1;
    node_type starting_node = 0;
    util::PriorityQueueParameters pq_params;
};

#ifdef COUNT_STATS
struct StatCounters {
    std::size_t pushed_nodes = 0;
    std::size_t ignored_nodes = 0;
    std::size_t extracted_nodes = 0;
    std::size_t failed_extracts = 0;
    std::size_t processed_nodes = 0;
    std::size_t idle = 0;
};

Settings settings;
Graph graph;
std::unique_ptr<StatCounters[]> stats;

inline void pushed_node(unsigned int id) { ++stats[id].pushed_nodes; }
inline void ignored_node(unsigned int id) { ++stats[id].ignored_nodes; }
inline void extracted_node(unsigned int id) { ++stats[id].extracted_nodes; }
inline void extract_failed(unsigned int id) { ++stats[id].failed_extracts; }
inline void processed_node(unsigned int id) { ++stats[id].processed_nodes; }
inline void started_idling(unsigned int id) { ++stats[id].idle; }
#else
inline void pushed_node(unsigned int) {}
inline void ignored_node(unsigned int) {}
inline void extracted_node(unsigned int) {}
inline void extract_failed(unsigned int) {}
inline void processed_node(unsigned int) {}
inline void started_idling(unsigned int) {}
#endif

// Each thread has a state, which is either working (0), check_idle (1) or idle
// (2) or woken up by another thread (3)
struct alignas(2 * L1_CACHE_LINESIZE) IdleState {
    std::atomic_uint state = 0;
};

std::unique_ptr<IdleState[]> idle_state;

alignas(2 * L1_CACHE_LINESIZE) std::atomic_size_t idle_counter = 0;

inline bool idle(unsigned int id) {
    idle_state[id].state.store(2, std::memory_order_relaxed);
    idle_counter.fetch_add(1, std::memory_order_relaxed);
    started_idling(id);
    while (true) {
        if (idle_state[id].state.load(std::memory_order_acquire) == 0) {
            // Someone woke this thread up
            return false;
        }
        if (idle_counter.load(std::memory_order_relaxed) ==
            2 * settings.num_threads) {
            // Everyone idles
            return true;
        }
        std::this_thread::yield();
    }
}

std::unique_ptr<Distance[]> shortest_distances;

static constexpr auto retries = 400;

// Returns true if this node inserted new nodes
bool process_node(PriorityQueue::Handle& handle,
                  PriorityQueue::value_type const& value, unsigned int id) {
    auto current_distance =
        shortest_distances[value.second].value.load(std::memory_order_relaxed);
    if (value.first > current_distance) {
        ignored_node(id);
        return false;
    }
    processed_node(id);
    bool inserted = false;
    for (node_type i = graph.nodes[value.second];
         i < graph.nodes[value.second + 1]; ++i) {
        node_type target = graph.edges[i].target;
        distance_type d = value.first + graph.edges[i].weight;
        distance_type old_d =
            shortest_distances[target].value.load(std::memory_order_relaxed);
        while (d < old_d) {
            if (shortest_distances[target].value.compare_exchange_weak(
                    old_d, d, std::memory_order_relaxed,
                    std::memory_order_relaxed)) {
                handle.push({d, target});
                pushed_node(id);
                inserted = true;
                break;
            }
        }
    }
    return inserted;
}

void main_loop(typename PriorityQueue::Handle& handle, unsigned int id) {
    auto extract_node = [&handle, id](auto& retval) {
        for (std::size_t i = 0; i < retries * settings.num_threads; ++i) {
            if (handle.try_extract_top(retval)) {
                return true;
            }
            extract_failed(id);
        }
        return false;
    };

#if defined PQ_MQ || defined PQ_MF
    auto partition_empty = [&handle, id]() {
        for (unsigned int i = 0; i < *settings.pq_params.c; ++i) {
            if (!handle.is_empty(id * (*settings.pq_params.c) + i)) {
                return false;
            }
        }
        return true;
    };
#endif

    PriorityQueue::value_type retval;

    while (true) {
        if (!extract_node(retval)) {
            // no item found, initiate idling
            idle_state[id].state.store(1, std::memory_order_relaxed);
            idle_counter.fetch_add(1, std::memory_order_release);
#if defined PQ_MQ || defined PQ_MF
            if (partition_empty()) {
                if (idle(id)) {
                    return;
                }
            } else {
                idle_state[id].state.store(0, std::memory_order_relaxed);
                idle_counter.fetch_sub(1, std::memory_order_release);
            }
            continue;
#else
            if (!handle.try_extract_top(retval)) {
                if (idle(id)) {
                    return;
                } else {
                    continue;
                }
            }
            idle_state[id].state.store(0, std::memory_order_relaxed);
            idle_counter.fetch_sub(1, std::memory_order_release);
#endif
        }
        extracted_node(id);
        if (process_node(handle, retval, id)) {
            if (idle_counter.load(std::memory_order_acquire) != 0) {
                for (std::size_t i = 0; i < settings.num_threads; ++i) {
                    if (i == id) {
                        continue;
                    }
                    unsigned int thread_idle_state =
                        idle_state[i].state.load(std::memory_order_relaxed);
                    while (true) {
                        while (thread_idle_state == 1) {
                            _mm_pause();
                            thread_idle_state = idle_state[i].state.load(
                                std::memory_order_relaxed);
                        }
                        if (thread_idle_state != 2) {
                            break;
                        }
                        if (idle_state[i].state.compare_exchange_strong(
                                thread_idle_state, 3,
                                std::memory_order_relaxed)) {
                            idle_counter.fetch_sub(2,
                                                   std::memory_order_relaxed);
                            idle_state[i].state.store(
                                0, std::memory_order_release);
                            break;
                        }
                    }
                }
            }
        }
    }
}

struct Task {
    static void run(thread_coordination::Context ctx, PriorityQueue& pq) {
        auto handle = pq.get_handle();

        if (ctx.is_main()) {
            std::clog << "Solving knapsack instance..." << std::flush;
            for (std::size_t i = 0; i + 1 < graph.nodes.size(); ++i) {
                shortest_distances[i].value =
                    std::numeric_limits<distance_type>::max();
            }
            shortest_distances[settings.starting_node].value.store(
                0, std::memory_order_relaxed);
            process_node(handle, {0, settings.starting_node}, ctx.get_id());
        }
        ctx.execute_synchronized_timed("main_loop", main_loop, handle,
                                       ctx.get_id());
        if (ctx.is_main()) {
            std::clog << "done\n";
        }
    }

    static threading::thread_config get_config(
        thread_coordination::Context const& ctx) {
        threading::thread_config config;
        config.cpu_set.reset();
        config.cpu_set.set(ctx.get_id());
        return config;
    }
};

void read_graph() {
    std::ifstream graph_stream{settings.graph_file};
    if (!graph_stream) {
        throw std::runtime_error{"Could not open graph file"};
    }
    std::vector<std::vector<Graph::Edge>> edges_per_node;
    node_type source;
    node_type target;
    distance_type weight;
    std::string problem;
    std::string first;
    while (graph_stream >> first) {
        if (first == "c") {
            graph_stream.ignore(std::numeric_limits<std::streamsize>::max(),
                                '\n');
        } else if (first == "p") {
            graph_stream >> problem;
            std::size_t num_nodes;
            std::size_t num_edges;
            graph_stream >> num_nodes >> num_edges;
            graph.nodes.resize(num_nodes + 1, 0);
            edges_per_node.resize(num_nodes);
            graph.edges.reserve(num_edges);
        } else if (first == "a") {
            graph_stream >> source >> target >> weight;
            edges_per_node[source - 1].push_back(
                Graph::Edge{target - 1, weight});
        } else {
            throw std::runtime_error{"Error reading file"};
        }
    }
    for (std::size_t i = 0; i < edges_per_node.size(); ++i) {
        graph.nodes[i + 1] =
            graph.nodes[i] + static_cast<node_type>(edges_per_node[i].size());
        std::copy(edges_per_node[i].begin(), edges_per_node[i].end(),
                  std::back_inserter(graph.edges));
    }
}

int main(int argc, char* argv[]) {
    std::clog << "Build configuration:\n\n";
#ifndef NDEBUG
    std::clog << "DEBUG build\n";
#endif
#ifdef COUNT_STATS
    std::clog << "Counting stats\n";
#endif
    std::clog << "L1 cache linesize (byte): " << L1_CACHE_LINESIZE << "\n";
    std::clog << '\n';

    std::clog << "Command line: ";
    std::copy(argv, argv + argc,
              std::ostream_iterator<char const*>(std::clog, " "));
    std::clog << "\n\n";

    cxxopts::Options options(
        "Shortest path benchmark",
        "This executable measures and records the performance of relaxed "
        "priority queues in the SSSP problem");
    // clang-format off
    options.add_options()
      ("j,threads", "Specify the number of threads "
       "(default: 4)", cxxopts::value<unsigned int>(), "NUMBER")
      ("n,start", "The starting node"
       "(default: 0)", cxxopts::value<node_type>(), "NUMBER")
      ("f,file", "The input graph", cxxopts::value<std::filesystem::path>(settings.graph_file)->default_value("graph.gr"), "PATH")
      ("s,seed", "Specify the initial seed"
       "(default: 0)", cxxopts::value<std::uint32_t>(), "NUMBER")
      ("o,output", "File to output the solution to (default: no output)", cxxopts::value<std::filesystem::path>(settings.output), "PATH")
      ("c,factor", "The number of queues when using multiqueue or multififo"
       "(default: 4)", cxxopts::value<std::size_t>(), "NUMBER")
      ("k,stickiness", "The stickiness when using multiqueue or multififo supporting stickiness"
       "(default: 8)", cxxopts::value<unsigned int>(), "NUMBER")
      ("h,help", "Print this help");
    // clang-format on

    try {
        auto result = options.parse(argc, argv);
        if (result.count("help")) {
            std::cerr << options.help() << std::endl;
            return 0;
        }
        if (result.count("start") > 0) {
            settings.starting_node = result["start"].as<node_type>();
        }
        if (result.count("threads") > 0) {
            settings.num_threads = result["threads"].as<unsigned int>();
        }
        if (result.count("factor") > 0) {
            settings.pq_params.c = result["factor"].as<std::size_t>();
        }
        if (result.count("stickiness") > 0) {
            settings.pq_params.stickiness =
                result["stickiness"].as<unsigned int>();
        }
        if (result.count("seed") > 0) {
            settings.seed = result["seed"].as<std::uint32_t>();
        }
    } catch (cxxopts::OptionParseException const& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

#ifndef NDEBUG
    std::clog << "Using debug build!\n\n";
#endif
    std::clog << "Settings: \n\t"
              << "Threads: " << settings.num_threads << "\n\t"
              << "Graph file: " << settings.graph_file.string() << "\n\t"
              << "Starting node: " << settings.starting_node << "\n\t"
              << "Seed: " << settings.seed;
    std::clog << "\n\n";

    xoroshiro256starstar rng;
    rng.seed(settings.seed);
    settings.pq_params.seed = rng();
    auto pq = util::create_pq<PriorityQueue>(1'000'000, settings.num_threads,
                                             settings.pq_params);

    std::clog << "Using priority queue: " << pq.description() << "\n\n";

    std::clog << "Reading graph..." << std::flush;
    try {
        read_graph();
    } catch (std::runtime_error const& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
    assert(graph.nodes.size() > 0);

    std::clog << "done\n";

    std::cout << "nodes: " << graph.nodes.size() - 1 << '\n';
    std::cout << "edges: " << graph.edges.size() << '\n';

    idle_counter = 0;
    idle_state = std::make_unique<IdleState[]>(settings.num_threads);
#ifdef COUNT_STATS
    stats = std::make_unique<StatCounters[]>(settings.num_threads);
#endif

    shortest_distances = std::make_unique<Distance[]>(graph.nodes.size() - 1);

    thread_coordination::ThreadCoordinator coordinator{settings.num_threads};
    coordinator.run_task<Task>(std::ref(pq));
    coordinator.join();
    auto duration = *coordinator.get_duration("main_loop");
    std::cout << "time: " << std::setprecision(3)
              << std::chrono::duration<double>(duration).count() << '\n';
#ifdef COUNT_STATS
    StatCounters total_stats = std::reduce(
        stats.get(), stats.get() + settings.num_threads, StatCounters{},
        [](StatCounters lhs, StatCounters const& rhs) {
            lhs.pushed_nodes += rhs.pushed_nodes;
            lhs.ignored_nodes += rhs.ignored_nodes;
            lhs.extracted_nodes += rhs.extracted_nodes;
            lhs.failed_extracts += rhs.failed_extracts;
            lhs.processed_nodes += rhs.processed_nodes;
            lhs.idle += rhs.idle;
            return lhs;
        });
    std::cout << "pushed nodes: " << total_stats.pushed_nodes << '\n';
    std::cout << "ignored nodes: " << total_stats.ignored_nodes << '\n';
    std::cout << "extracted nodes: " << total_stats.extracted_nodes << '\n';
    std::cout << "failed extracts: " << total_stats.failed_extracts << '\n';
    std::cout << "processed nodes: " << total_stats.processed_nodes << '\n';
    std::cout << "idle: " << total_stats.idle << '\n';
#endif
    if (!settings.output.empty()) {
        std::clog << "Writing output..." << std::flush;
        std::ofstream out_stream{settings.output};
        if (!out_stream) {
            std::cerr << "Could not open output file\n";
            return 1;
        }
        out_stream << "node dist\n";
        for (std::size_t i = 0; i + 1 < graph.nodes.size(); ++i) {
            out_stream << i << ' '
                       << shortest_distances[i].value.load(
                              std::memory_order_relaxed)
                       << '\n';
        }
        std::clog << "done\n";
    }
}
