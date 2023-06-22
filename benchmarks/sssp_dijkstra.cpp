#include "graph.hpp"
#include "priority_queue_factory.hpp"
#include "termination_detection.hpp"
#include "thread_coordination.hpp"

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

using PriorityQueue = DefaultMinPriorityQueue;

using distance_type = PriorityQueue::key_type;
using node_type = PriorityQueue::mapped_type;
using handle_type = PriorityQueue::Handle;

using thread_coordination::timepoint_type;

struct Settings {
    int num_threads = 4;
    std::filesystem::path graph_file;
    std::filesystem::path solution_file;
    std::filesystem::path distance_file;
    int seed = 1;
};

struct ThreadData {
    handle_type pq_handle;
    long long pushed_nodes{0};
    long long ignored_nodes{0};
    long long processed_nodes{0};
};

struct SharedData {
    template <typename T>
    struct alignas(L1_CACHE_LINESIZE) PaddedAtomic {
        std::atomic<T> value;
    };
    std::vector<PaddedAtomic<distance_type>> shortest_distances;
    std::atomic<timepoint_type> start_time{timepoint_type::max()};
    std::atomic<timepoint_type> end_time{timepoint_type::min()};
    std::atomic_llong pushed_nodes{0};
    std::atomic_llong ignored_nodes{0};
    std::atomic_llong processed_nodes{0};
    termination_detection::Data termination_detection_data;

    SharedData(std::size_t num_nodes) : shortest_distances(num_nodes) {
        for (auto& i : shortest_distances) {
            i.value.store(std::numeric_limits<distance_type>::max(), std::memory_order_relaxed);
        }
    }

    void update_work_time(std::pair<timepoint_type, timepoint_type> const& work_time) {
        auto old_start = start_time.load(std::memory_order_relaxed);
        while (work_time.first < old_start &&
               !start_time.compare_exchange_weak(old_start, work_time.first, std::memory_order_relaxed)) {
        }
        auto old_end = start_time.load(std::memory_order_relaxed);
        while (work_time.second > old_end &&
               !end_time.compare_exchange_weak(old_end, work_time.second, std::memory_order_relaxed)) {
        }
    }

    void update_counters(ThreadData const& thread_data) {
        pushed_nodes.fetch_add(thread_data.pushed_nodes, std::memory_order_relaxed);
        ignored_nodes.fetch_add(thread_data.ignored_nodes, std::memory_order_relaxed);
        processed_nodes.fetch_add(thread_data.processed_nodes, std::memory_order_relaxed);
    }
};

void process_node(PriorityQueue::value_type const& value, SharedData& shared_data, ThreadData& data,
                  Graph const& graph) {
    auto current_distance = shared_data.shortest_distances[value.second].value.load(std::memory_order_relaxed);
    if (value.first > current_distance) {
        ++data.ignored_nodes;
        return;
    }
    ++data.processed_nodes;
    for (node_type i = graph.nodes[value.second]; i < graph.nodes[value.second + 1]; ++i) {
        node_type target = graph.edges[i].target;
        distance_type d = value.first + graph.edges[i].weight;
        distance_type old_d = shared_data.shortest_distances[target].value.load(std::memory_order_relaxed);
        while (d < old_d) {
            if (shared_data.shortest_distances[target].value.compare_exchange_weak(old_d, d, std::memory_order_seq_cst,
                                                                                   std::memory_order_relaxed)) {
                data.pq_handle.push({d, target});
                ++data.pushed_nodes;
                break;
            }
        }
    }
}

void main_loop(int num_threads, SharedData& shared_data, ThreadData& data, Graph const& graph) {
    PriorityQueue::value_type retval;
    while (termination_detection::try_do(num_threads, shared_data.termination_detection_data,
                                         [&]() { return data.pq_handle.try_pop(retval); })) {
        process_node(retval, shared_data, data, graph);
    }
}

void benchmark_thread(thread_coordination::Context ctx, PriorityQueue& pq, SharedData& shared_data,
                      Graph const& graph) {
    ThreadData data{pq.get_handle(ctx.get_id())};
    if (ctx.get_id() == 0) {
        shared_data.shortest_distances[0].value = 0;
        data.pq_handle.push({0, 0});
        ++data.pushed_nodes;
    }
    auto work_time = ctx.execute_synchronized([&]() { main_loop(ctx.get_num_threads(), shared_data, data, graph); });
    shared_data.update_work_time(work_time);
    shared_data.update_counters(data);
}

void print_settings(Settings const& settings) {
    std::clog << "Threads: " << settings.num_threads << '\n'
              << "Graph: " << settings.graph_file.string() << '\n'
              << "Seed: " << settings.seed;
    std::clog << "\n\n";
}

void print_shared_data(Settings const& settings, SharedData const& shared_data, Graph const& graph, bool valid) {
    auto time = std::chrono::duration<double>(shared_data.end_time.load() - shared_data.start_time.load()).count();
    std::clog << "Time (s): " << std::setprecision(3) << time << '\n';
    std::clog << "Processed nodes: " << shared_data.processed_nodes.load() << '\n';
    std::clog << "Ignored nodes: " << shared_data.ignored_nodes.load() << '\n';
    std::clog << "Total nodes: " << shared_data.processed_nodes.load() + shared_data.ignored_nodes.load() << '\n';

    std::cout << "graph,nodes,edges,threads,seed,work_time,processed_nodes,ignored_nodes,valid\n";
    std::cout << settings.graph_file.string() << ',' << graph.num_nodes() << ',' << graph.num_edges() << ','
              << settings.num_threads << ',' << settings.seed << ',' << time << ','
              << shared_data.processed_nodes.load() << ',' << shared_data.ignored_nodes << ',' << valid << '\n';
}

bool verify_distances(std::filesystem::path const& file, SharedData const& shared_data) {
    std::ifstream in{file};
    if (!in) {
        throw std::runtime_error("Failed to open file: " + file.string());
    }
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream iss{line};
        node_type node = 0;
        distance_type distance = 0;
        iss >> node >> distance;

        if (shared_data.shortest_distances[node].value != distance) {
            return false;
        }
    }
    return true;
}

void write_distances(SharedData const& shared_data, std::filesystem::path const& file) {
    std::ofstream out{file};
    if (!out) {
        throw std::runtime_error("Failed to open file: " + file.string());
    }
    for (std::size_t i = 0; i < shared_data.shortest_distances.size(); ++i) {
        out << i << ' ' << shared_data.shortest_distances[i].value << '\n';
    }
}

bool run_benchmark(Settings const& settings, PriorityQueueConfig const& pq_config) {
    std::clog << "Reading graph..." << std::flush;
    Graph graph(0, 0);
    try {
        graph.from_file(settings.graph_file);
    } catch (std::runtime_error const& e) {
        std::clog << "failed: " << e.what() << std::endl;
        return false;
    }
    std::clog << "done\n";

    auto pq = create_pq<PriorityQueue>(settings.num_threads, graph.num_nodes(), pq_config);
    SharedData shared_data{graph.num_nodes()};
    std::clog << "Computing shortest paths..." << std::flush;
    thread_coordination::TaskHandle task_handle{settings.num_threads, benchmark_thread, std::ref(pq),
                                                std::ref(shared_data), graph};
    task_handle.wait();

    std::clog << "done" << std::endl;

    bool success = true;
    if (!settings.distance_file.empty()) {
        std::clog << "Writing distances..." << std::flush;
        try {
            write_distances(shared_data, settings.distance_file);
            std::clog << "done" << std::endl;
        } catch (std::runtime_error const& e) {
            std::clog << "failed: " << e.what() << std::endl;
            success = false;
        }
    }
    std::clog << "Verifying..." << std::flush;
    bool valid = true;
    if (shared_data.processed_nodes.load() + shared_data.ignored_nodes.load() != shared_data.pushed_nodes.load()) {
        std::clog << "failed: Not all nodes were popped" << std::endl;
        success = false;
        valid = false;
    } else {
        try {
            if (verify_distances(settings.solution_file, shared_data)) {
                std::clog << "done" << std::endl;
            } else {
                std::clog << "failed: Distance mismatch" << std::endl;
                success = false;
                valid = false;
            }
        } catch (std::runtime_error const& e) {
            std::clog << "failed: " << e.what() << std::endl;
            valid = false;
            success = false;
        }
    }
    std::clog << '\n';
    print_shared_data(settings, shared_data, graph, valid);
    return success;
}

void print_header() {
    std::clog << "Built on " << __DATE__ << ' ' << __TIME__ << " with:\n";
#ifdef NDEBUG
    std::clog << "  Release build\n";
#else
    std::clog << "  Debug build\n";
#endif
#ifdef __clang__
    std::clog << "  Clang " << __clang_version__ << '\n';
#elif defined(__GNUC__)
    std::clog << "  GCC " << __VERSION__ << '\n';
#else
    std::clog << "  Unknown compiler\n";
#endif
    std::clog << "  Priority queue: " << pq_name << '\n';
}

int main(int argc, char* argv[]) {
    print_header();
    std::clog << "\nCommand line:";
    for (int i = 0; i < argc; ++i) {
        std::clog << ' ' << argv[i];
    }
    std::clog << '\n' << '\n';

    Settings settings{};
    cxxopts::Options options(argv[0]);
    // clang-format off
    options.add_options()
      ("j,threads", "The number of threads", cxxopts::value<int>(settings.num_threads), "NUMBER")
      ("file", "The input graph", cxxopts::value<std::filesystem::path>(settings.graph_file), "PATH")
      ("solution", "The reference solution", cxxopts::value<std::filesystem::path>(settings.solution_file), "PATH")
      ("o,distance-file", "Path to write the distances to", cxxopts::value<std::filesystem::path>(settings.distance_file), "PATH")
      ("h,help", "Print this help");
    // clang-format on
    add_pq_options(options);
    options.parse_positional({"file", "solution"});

    PriorityQueueConfig pq_config;
    {
        cxxopts::ParseResult result;
        try {
            result = options.parse(argc, argv);
        } catch (cxxopts::OptionParseException const& e) {
            std::cerr << "Error parsing arguments: " << e.what() << '\n';
            std::cerr << options.help() << std::endl;
            return 1;
        }
        if (result.count("help") > 0) {
            std::cerr << options.help() << std::endl;
            return 0;
        }
        pq_config = get_pq_options(result);
    }
    if (settings.graph_file.empty()) {
        std::cerr << "Error: No graph file specified" << std::endl;
        std::cerr << options.help() << std::endl;
        return 1;
    }
    if (settings.solution_file.empty()) {
        std::cerr << "Error: No solution file specified" << std::endl;
        std::cerr << options.help() << std::endl;
        return 1;
    }

    print_settings(settings);

    bool success = run_benchmark(settings, pq_config);

    return success ? 0 : 1;
}
