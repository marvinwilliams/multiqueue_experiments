#include "priority_queue_factory.hpp"
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

struct Graph {
    struct Edge {
        node_type target;
        distance_type weight;
    };
    std::vector<node_type> nodes;
    std::vector<Edge> edges;
};

struct ThreadData {
    handle_type pq_handle;
    long long pushed_nodes{0};
    long long ignored_nodes{0};
    long long popped_nodes{0};
    long long processed_nodes{0};
    long long idling{0};
};

struct SharedData {
    template <typename T>
    struct alignas(L1_CACHE_LINESIZE) PaddedAtomic {
        std::atomic<T> value;
    };
    alignas(L1_CACHE_LINESIZE) std::atomic_int no_work_count{0};
    int idle_count{0};
    std::mutex idle_mutex;
    std::condition_variable idle_cv;
    std::vector<PaddedAtomic<distance_type>> shortest_distances;
    std::atomic<timepoint_type> start_time{timepoint_type::max()};
    std::atomic<timepoint_type> end_time{timepoint_type::min()};
    std::atomic_llong pushed_nodes{0};
    std::atomic_llong ignored_nodes{0};
    std::atomic_llong popped_nodes{0};
    std::atomic_llong processed_nodes{0};

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
        popped_nodes.fetch_add(thread_data.popped_nodes, std::memory_order_relaxed);
        processed_nodes.fetch_add(thread_data.processed_nodes, std::memory_order_relaxed);
    }
};

Graph read_graph(std::filesystem::path const& graph_file) {
    int fd = open(graph_file.c_str(), O_RDONLY);
    if (fd == -1) {
        throw std::runtime_error{"Could not open file"};
    }

    struct stat sb {};

    if (fstat(fd, &sb) == -1) {
        close(fd);
        throw std::runtime_error{"Could not get file size"};
    }

    auto* addr = mmap(nullptr, static_cast<std::size_t>(sb.st_size), PROT_READ, MAP_PRIVATE, fd, 0);
    if (addr == MAP_FAILED) {
        close(fd);
        throw std::runtime_error{"mmap failed"};
    }
    close(fd);
    madvise(addr, static_cast<std::size_t>(sb.st_size), MADV_SEQUENTIAL);

    std::vector<std::pair<std::size_t, Graph::Edge>> edge_list;
    std::size_t num_nodes = 0;
    std::size_t num_edges = 0;
    const auto* it = static_cast<char const*>(addr);
    const auto* end = it + sb.st_size;
    while (*it == 'c') {
        ++it;
        while (*it++ != '\n') {
        }
    }
    Graph graph;
    if (*it == 'p') {
        while (std::isspace(*++it) != 0) {
        }
        while (std::isspace(*++it) == 0) {
        }
        while (std::isspace(*++it) != 0) {
        }
        auto res = std::from_chars(it, end, num_nodes);
        assert(res.ec == std::errc{});
        it = res.ptr;
        while (std::isspace(*++it) != 0) {
        }
        res = std::from_chars(it, end, num_edges);
        assert(res.ec == std::errc{});
        it = res.ptr;
        graph.nodes.resize(num_nodes + 1);
        graph.edges.resize(num_edges);
        edge_list.reserve(num_edges);
        while (std::isspace(*it) != 0) {
            ++it;
        }
    }
    while (*it == 'c') {
        ++it;
        while (*it++ != '\n') {
        }
    }
    std::size_t source = 0;
    std::size_t target = 0;
    unsigned int weight = 0;
    for (std::size_t i = 0; i != num_edges; ++i) {
        while (std::isspace(*it) != 0) {
            ++it;
        }
        assert(*it == 'a' && std::isspace(*(it + 1)));
        it += 2;
        auto res = std::from_chars(it, end, source);
        assert(res.ec == std::errc{});
        it = res.ptr;
        while (std::isspace(*it) != 0) {
            ++it;
        }
        res = std::from_chars(it, end, target);
        assert(res.ec == std::errc{});
        it = res.ptr;
        while (std::isspace(*it) != 0) {
            ++it;
        }
        res = std::from_chars(it, end, weight);
        assert(res.ec == std::errc{});
        it = res.ptr;
        // 1-based
        edge_list.push_back({source - 1, {target - 1, weight}});
    }
    munmap(addr, static_cast<std::size_t>(sb.st_size));
    std::sort(edge_list.begin(), edge_list.end(), [](auto& lhs, auto& rhs) { return lhs.first < rhs.first; });
    graph.nodes[0] = 0;
    for (std::size_t i = 1; i < num_nodes + 1; ++i) {
        graph.nodes[i] = graph.nodes[i - 1];
        while (graph.nodes[i] < num_edges && edge_list[graph.nodes[i]].first == i - 1) {
            ++graph.nodes[i];
        }
    }
    std::transform(edge_list.begin(), edge_list.end(), graph.edges.begin(), [](auto const& e) { return e.second; });
    return graph;
}

// Returns true if this node inserted new nodes
bool process_node(PriorityQueue::value_type const& value, SharedData& shared_data, ThreadData& data,
                  Graph const& graph) {
    auto current_distance = shared_data.shortest_distances[value.second].value.load(std::memory_order_relaxed);
    if (value.first > current_distance) {
        ++data.ignored_nodes;
        return false;
    }
    ++data.processed_nodes;
    bool inserted = false;
    for (node_type i = graph.nodes[value.second]; i < graph.nodes[value.second + 1]; ++i) {
        node_type target = graph.edges[i].target;
        distance_type d = value.first + graph.edges[i].weight;
        distance_type old_d = shared_data.shortest_distances[target].value.load(std::memory_order_relaxed);
        while (d < old_d) {
            if (shared_data.shortest_distances[target].value.compare_exchange_weak(old_d, d, std::memory_order_seq_cst,
                                                                                   std::memory_order_relaxed)) {
                data.pq_handle.push({d, target});
                ++data.pushed_nodes;
                inserted = true;
                break;
            }
        }
    }
    return inserted;
}

void benchmark_thread(thread_coordination::Context ctx, PriorityQueue& pq, SharedData& shared_data,
                      Graph const& graph) {
    ThreadData data{pq.get_handle(ctx.get_id())};
    if (ctx.get_id() == 0) {
        shared_data.shortest_distances[0].value = 0;
        data.pq_handle.push({0, 0});
        ++data.pushed_nodes;
    }
    auto work_time = ctx.execute_synchronized([&]() {
        while (true) {
            PriorityQueue::value_type retval;
            if (!data.pq_handle.try_pop(retval)) {
                ctx.write(std::cerr) << "No work found\n";
                auto num_no_work = shared_data.no_work_count.fetch_add(1, std::memory_order_relaxed) + 1;
                while (num_no_work < ctx.get_num_threads()) {
                    if (data.pq_handle.try_pop(retval)) {
                        shared_data.no_work_count.fetch_sub(1, std::memory_order_relaxed);
                        ctx.write(std::cerr) << "Found work\n";
                        break;
                    }
                    num_no_work = shared_data.no_work_count.load(std::memory_order_relaxed);
                }
                if (num_no_work >= ctx.get_num_threads()) {
                    ctx.write(std::cerr) << "I think all threads finished, final check\n";
                    if (data.pq_handle.try_pop(retval)) {
                        shared_data.no_work_count.fetch_sub(1, std::memory_order_relaxed);
                        ctx.write(std::cerr) << "Found work\n";
                    } else {
                        ctx.write(std::cerr) << "Idling now\n";
                        auto l = std::unique_lock(shared_data.idle_mutex);
                        if (++shared_data.idle_count == ctx.get_num_threads()) {
                            shared_data.idle_cv.notify_all();
                            ctx.write(std::cerr) << "Last worker thread finished\n";
                            return;
                        }
                        shared_data.idle_cv.wait(l, [&]() {
                            return shared_data.idle_count == ctx.get_num_threads() ||
                                shared_data.no_work_count.load(std::memory_order_relaxed) < ctx.get_num_threads();
                        });
                        if (shared_data.idle_count == ctx.get_num_threads()) {
                            ctx.write(std::cerr) << "Worker thread finished\n";
                            return;
                        }
                        --shared_data.idle_count;
                        l.unlock();
                        shared_data.no_work_count.fetch_sub(1, std::memory_order_relaxed);
                        ctx.write(std::cerr) << "Others are working again\n";
                        continue;
                    }
                }
            }
            ++data.popped_nodes;
            process_node(retval, shared_data, data, graph);
        }
    });
    shared_data.update_work_time(work_time);
    shared_data.update_counters(data);
}

void run_benchmark(Settings const& settings, PriorityQueueConfig const& pq_config, SharedData& shared_data,
                   Graph const& graph) {
    auto pq = create_pq<PriorityQueue>(settings.num_threads, graph.nodes.size() - 1, pq_config);

    std::clog << "Computing shortest paths..." << std::flush;

    thread_coordination::TaskHandle task_handle{settings.num_threads, benchmark_thread, std::ref(pq),
                                                std::ref(shared_data), graph};
    task_handle.wait();

    std::clog << "done" << std::endl;
}

bool verify_shared_data(Settings const& settings, SharedData const& shared_data) {
    std::ifstream in{settings.solution_file};
    if (!in) {
        throw std::runtime_error("Failed to open file: " + settings.solution_file.string());
    }
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream iss{line};
        node_type node;
        distance_type distance;
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

void print_settings(Settings const& settings) {
    std::clog << "Threads: " << settings.num_threads << '\n'
              << "Graph: " << settings.graph_file.string() << '\n'
              << "Seed: " << settings.seed;
    std::clog << "\n\n";
}

void print_shared_data(Settings const& settings, SharedData const& shared_data, Graph const& graph) {
    std::clog << "Time (s): " << std::setprecision(3)
              << std::chrono::duration<double>(shared_data.end_time.load() - shared_data.start_time.load()).count()
              << '\n';
    std::cout << "Pushed nodes: " << shared_data.pushed_nodes << '\n';
    std::cout << "Ignored nodes: " << shared_data.ignored_nodes << '\n';
    std::cout << "Popped nodes: " << shared_data.popped_nodes << '\n';
    std::cout << "Processed nodes: " << shared_data.processed_nodes << '\n';

    std::cout << "threads,graph,nodes,edges,seed,work_time,pushed_nodes,ignored_nodes,popped_nodes,"
                 "processed_nodes\n";
    std::cout << settings.num_threads << ',' << settings.graph_file.string() << ',' << graph.nodes.size() - 1 << ','
              << graph.edges.size() << ',' << settings.seed << ','
              << std::chrono::duration<double>(shared_data.end_time.load() - shared_data.start_time.load()).count()
              << ',' << shared_data.pushed_nodes << ',' << shared_data.ignored_nodes << ',' << shared_data.popped_nodes
              << ',' << shared_data.processed_nodes << '\n';
}

int main(int argc, char* argv[]) {
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
    std::clog << '\n';

    std::clog << "Command line:";
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
        cxxopts::ParseResult args;
        try {
            args = options.parse(argc, argv);
        } catch (cxxopts::OptionParseException const& e) {
            std::cerr << "Error parsing arguments: " << e.what() << '\n';
            std::cerr << options.help() << std::endl;
            return 1;
        }
        if (args.count("help") > 0) {
            std::cerr << options.help() << std::endl;
            return 0;
        }
        pq_config = get_pq_options(args);
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

    std::clog << "Reading graph..." << std::flush;
    Graph graph;
    try {
        graph = read_graph(settings.graph_file);
    } catch (std::runtime_error const& e) {
        std::cerr << "\nError reading graph: " << e.what() << std::endl;
        return 1;
    }
    assert(!graph.nodes.empty());
    std::clog << "done\n";
    std::clog << '\n';

    SharedData shared_data{graph.nodes.size() - 1};
    run_benchmark(settings, pq_config, shared_data, graph);

    if (!settings.distance_file.empty()) {
        std::clog << "Writing distances..." << std::flush;
        try {
            write_distances(shared_data, settings.distance_file);
            std::clog << "done" << std::endl;
        } catch (std::runtime_error const& e) {
            std::clog << "failed: " << e.what() << std::endl;
        }
    }
    std::clog << "Verifying..." << std::flush;
    try {
        if (!verify_shared_data(settings, shared_data)) {
            std::clog << "failed: Distance mismatch" << std::endl;
            return 1;
        }
        std::clog << "done" << std::endl;
    } catch (std::runtime_error const& e) {
        std::clog << "failed: " << e.what() << std::endl;
        return 1;
    }

    std::clog << '\n';
    print_shared_data(settings, shared_data, graph);

    return 0;
}
