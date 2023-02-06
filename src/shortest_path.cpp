#include "thread_coordination.hpp"
#include "utils/priority_queue_factory.hpp"
#include "utils/thread_coordination.hpp"

#ifdef PQ_MQ
#include "multiqueue/config.hpp"
#endif

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

using PriorityQueue = typename util::PriorityQueueFactory<distance_type, node_type>::type;

struct Graph {
    struct Edge {
        node_type target;
        distance_type weight;
    };
    std::vector<node_type> nodes;
    std::vector<Edge> edges;
};

struct ThreadData {
    typename PriorityQueue::handle_type pq_handle;
#ifdef EXP_COUNT_STATS
    unsigned long long pushed_nodes{0};
    unsigned long long ignored_nodes{0};
    unsigned long long extracted_nodes{0};
    unsigned long long failed_extracts{0};
    unsigned long long processed_nodes{0};
#endif

    explicit ThreadData(typename PriorityQueue::handle_type h) : pq_handle(std::move(h)){};
};

#ifdef MULTIQUEUE_COUNT_STATS
#define INC_COUNTER(counter) ++data.counter
#else
#define INC_COUNTER(counter) (void)0
#endif

struct Result {
    std::chrono::milliseconds duration;
#ifdef EXP_COUNT_STATS
    alignas(L1_CACHE_LINESIZE) std::atomic_ullong pushed_nodes{0};
    alignas(L1_CACHE_LINESIZE) std::atomic_ullong ignored_nodes{0};
    alignas(L1_CACHE_LINESIZE) std::atomic_ullong popped_nodes{0};
    alignas(L1_CACHE_LINESIZE) std::atomic_ullong failed_pops{0};
    alignas(L1_CACHE_LINESIZE) std::atomic_ullong processed_nodes{0};
#endif
};

Graph graph;
Result result;

void read_graph(std::filesystem::path const& graph_file) {
    int fd = open(graph_file.c_str(), O_RDONLY);
    if (fd == -1) {
        throw std::runtime_error{"Could not open file"};
    }

    struct stat sb;

    if (fstat(fd, &sb) == -1) {
        close(fd);
        throw std::runtime_error{"Could not get file size"};
    }

    auto addr = mmap(NULL, static_cast<std::size_t>(sb.st_size), PROT_READ, MAP_PRIVATE, fd, 0);
    if (addr == MAP_FAILED) {
        close(fd);
        throw std::runtime_error{"mmap failed"};
    }
    close(fd);
    madvise(addr, static_cast<std::size_t>(sb.st_size), MADV_SEQUENTIAL);

    std::vector<std::pair<std::size_t, Graph::Edge>> edge_list;
    std::size_t num_nodes = 0;
    std::size_t num_edges = 0;
    auto it = static_cast<char const*>(addr);
    auto end = it + sb.st_size;
    while (*it == 'c') {
        ++it;
        while (*it++ != '\n') {
        }
    }
    if (*it == 'p') {
        while (std::isspace(*++it)) {
        }
        while (!std::isspace(*++it)) {
        }
        while (std::isspace(*++it)) {
        }
        auto res = std::from_chars(it, end, num_nodes);
        assert(res.ec == std::errc{});
        it = res.ptr;
        while (std::isspace(*++it)) {
        }
        res = std::from_chars(it, end, num_edges);
        assert(res.ec == std::errc{});
        it = res.ptr;
        graph.nodes.resize(num_nodes + 1);
        graph.edges.resize(num_edges);
        edge_list.reserve(num_edges);
        while (std::isspace(*it)) {
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
        while (std::isspace(*it)) {
            ++it;
        }
        assert(*it == 'a' && std::isspace(*(it + 1)));
        it += 2;
        auto res = std::from_chars(it, end, source);
        assert(res.ec == std::errc{});
        it = res.ptr;
        while (std::isspace(*it)) {
            ++it;
        }
        res = std::from_chars(it, end, target);
        assert(res.ec == std::errc{});
        it = res.ptr;
        while (std::isspace(*it)) {
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
}

// Each thread has a state, which is either working (0), check_idle (1) or idle
// (2) or woken up by another thread (3)
struct alignas(2 * L1_CACHE_LINESIZE) IdleState {
    std::atomic_uint state = 0;
};

alignas(L1_CACHE_LINESIZE) std::atomic_uint idle_counter{0};
std::vector<IdleState> idle_state;

bool idle(unsigned int num_threads, unsigned int id) {
    idle_state[id].state.store(2, std::memory_order_relaxed);
    idle_counter.fetch_add(1, std::memory_order_relaxed);
    while (true) {
        if (idle_state[id].state.load(std::memory_order_acquire) == 0) {
            // Someone woke this thread up
            return false;
        }
        if (idle_counter.load(std::memory_order_relaxed) == 2 * num_threads) {
            // Everyone idles
            return true;
        }
        _mm_pause();
    }
}

struct alignas(2 * L1_CACHE_LINESIZE) Distance {
    std::atomic<distance_type> value;
};

std::vector<Distance> shortest_distances;

// Returns true if this node inserted new nodes
bool process_node(ThreadData& data, PriorityQueue::value_type const& value) {
    auto current_distance = shortest_distances[value.second].value.load(std::memory_order_relaxed);
    if (value.first > current_distance) {
        INC_COUNTER(ignored_nodes);
        return false;
    }
    INC_COUNTER(processed_nodes);
    bool inserted = false;
    for (node_type i = graph.nodes[value.second]; i < graph.nodes[value.second + 1]; ++i) {
        node_type target = graph.edges[i].target;
        distance_type d = value.first + graph.edges[i].weight;
        distance_type old_d = shortest_distances[target].value.load(std::memory_order_relaxed);
        while (d < old_d) {
            if (shortest_distances[target].value.compare_exchange_weak(old_d, d, std::memory_order_seq_cst,
                                                                       std::memory_order_relaxed)) {
                data.pq_handle.push({d, target});
                INC_COUNTER(pushed_nodes);
                inserted = true;
                break;
            }
        }
    }
    return inserted;
}

void main_loop(thread_coordination::Context ctx, ThreadData& data) {
    auto search_element = [&data](auto& retval) {
#if defined PQ_MQ || defined PQ_MF
        for (auto i = from; i != to; ++i) {
            if (data.pq_handle.try_pop_from(i, retval)) {
                return true;
            }
        }
        return false;
#else
        return handle.try_pop(retval);
#endif
    };

    PriorityQueue::value_type retval;

    while (true) {
        if (!handle.try_pop(retval)) {
            INC_COUNTER(pop_failed);
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
                for (std::size_t i = 0; i < num_threads; ++i) {
                    if (i == id) {
                        continue;
                    }
                    while (true) {
                        auto thread_idle_state = idle_state[i].state.load(std::memory_order_relaxed);
                        while (thread_idle_state == 1) {
                            _mm_pause();
                            thread_idle_state = idle_state[i].state.load(std::memory_order_relaxed);
                        }
                        if (thread_idle_state != 2) {
                            break;
                        }
                        if (idle_state[i].state.compare_exchange_strong(thread_idle_state, 3,
                                                                        std::memory_order_relaxed)) {
                            idle_counter.fetch_sub(2, std::memory_order_relaxed);
                            idle_state[i].state.store(0, std::memory_order_release);
                            break;
                        }
                    }
                }
            }
        }
    }
}

void run_benchmark(PriorityQueue& pq) {
    shortest_distances[settings.starting_node].value = 0;
    for (std::size_t i = 1; i < graph.nodes.size() - 1; ++i) {
        shortest_distances[i].value = std::numeric_limits<distance_type>::max();
    }

    std::vector std::clog << "\nComputing shortest paths..." << std::flush;

    thread_coordination::TaskHandletask_handle{settings.num_threads,
                                               [](thread_coordination::Context ctx, PriorityQueue& pq) {
                                                   ThreadData data{pq.get_handle(ctx.get_id())};
                                                   ctx.execute_synchronized_timed(duration, main_loop, data, ctx);
                                               },
                                               std::ref(pq)};
    task_handle.wait();
    std::clog << "done\n";
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
    std::copy(argv, argv + argc, std::ostream_iterator<char const*>(std::clog, " "));
    std::clog << "\n\n";

    std::filesystem::path graph_file;
    std::filesystem::path solution_file;
    std::filesystem::path out_file;
    unsigned int num_threads = 4;
    node_type starting_node = 0;

#ifdef PQ_MQ
    multiqueue::Config mq_config;
#endif

    cxxopts::Options options("Shortest path benchmark",
                             "This executable measures and records the performance of relaxed "
                             "priority queues in the SSSP problem");
    // clang-format off
    options.add_options()
      ("j,threads", "Specify the number of threads "
       "(default: 4)", cxxopts::value<unsigned int>(num_threads), "NUMBER")
      ("n,start", "The starting node"
       "(default: 0)", cxxopts::value<node_type>(starting_node), "NUMBER")
      ("f,file", "The input graph", cxxopts::value<std::filesystem::path>(graph_file)->default_value("graph.gr"), "PATH")
      ("v,solution", "The reference solution", cxxopts::value<std::filesystem::path>(solution_file)->default_value("solution.txt"), "PATH")
      ("o,output", "The output", cxxopts::value<std::filesystem::path>(out_file)->default_value("solution.txt"), "PATH")
#ifdef PQ_MQ
      ("s,seed", "Specify the initial seed"
       "(default: 0)", cxxopts::value<std::uint32_t>(mq_config.seed), "NUMBER")
      ("c,factor", "The number of queues when using multiqueue or multififo"
       "(default: 4)", cxxopts::value<std::size_t>(mq_config.c), "NUMBER")
      ("k,stickiness", "The stickiness when using multiqueue or multififo supporting stickiness"
       "(default: 8)", cxxopts::value<unsigned int>(mq_config.stickiness), "NUMBER")
#endif
      ("h,help", "Print this help");
    // clang-format on

    try {
        auto result = options.parse(argc, argv);
        if (result.count("help")) {
            std::cerr << options.help() << std::endl;
            return 0;
        }
    } catch (cxxopts::OptionParseException const& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    std::clog << "Settings: \n\t"
              << "Threads: " << num_threads << "\n\t"
              << "Graph file: " << graph_file.string() << "\n\t"
              << "Solution file: " << (solution_file ? solution_file.string() : "Not supplied") << "\n\t"
              << "Starting node: " << starting_node << "\n\t"
              << "Seed: " << settings.seed;
    std::clog << "\n\n";

#ifdef PQ_MQ
    auto pq = PriorityQueue(num_threads, mq_config);
#else
    auto pq = PriorityQueue(num_threads);
#endif

    std::clog << "Data structure: ";
    util::describe::describe(std::clog, pq);
    std::clog << '\n';

    std::clog << "Reading graph..." << std::flush;
    try {
        read_graph();
    } catch (std::runtime_error const& e) {
        std::cerr << '\n' << e.what() << '\n';
        return 1;
    }
    assert(graph.nodes.size() > 0);

    std::clog << "done\n";
    std::clog << "nodes: " << graph.nodes.size() - 1 << '\n';
    std::clog << "edges: " << graph.edges.size() << '\n';

    idle_state = std::vector<IdleState>(settings.num_threads);
    shortest_distances = std::vector<Distance>(graph.nodes.size() - 1);

    run_benchmark(pq);

    std::clog << "time: " << std::setprecision(3) << std::chrono::duration<double>(duration).count() << '\n';
#ifdef EXP_COUNT_STATS
    std::cout << "pushed nodes: " << result.pushed_nodes << '\n';
    std::cout << "ignored nodes: " << result.ignored_nodes << '\n';
    std::cout << "extracted nodes: " << result.popped_nodes << '\n';
    std::cout << "failed extracts: " << result.failed_pops << '\n';
    std::cout << "processed nodes: " << result.processed_nodes << '\n';
#endif
    if (!output.empty()) {
        std::ofstream out_stream{settings.output};
        if (out_stream) {
            std::clog << "\nWriting output..." << std::flush;
            out_stream << "node dist\n";
            for (std::size_t i = 0; i + 1 < graph.nodes.size(); ++i) {
                out_stream << i << ' ' << shortest_distances[i].value.load() << '\n';
            }
            std::clog << "done\n";
        } else {
            std::cerr << "\nCould not open output file\n";
        }
    }
    std::cout << num_threads << ',' << seed << ',' << std::fixed << std::setprecision(3)
              << std::chrono::duration<double>(result.duration).count();
#ifdef MULTIQUEUE_COUNT_STATS
    std::cout << ',' << result.pushed_nodes << ',' << result.igored_nodes << ',' << result.popped_nodes << ','
              << result.failed_pops << ',' << result.processed_nodes;
#else
    std::cout << ",n/a,n/a,n/a,n/a,n/a";
#endif
    std::cout << std::endl;
}
