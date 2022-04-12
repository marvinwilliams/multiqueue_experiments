#include "utils/priority_queue_factory.hpp"
#include "utils/thread_coordination.hpp"
#include "utils/threading.hpp"
#include "utils/xoroshiro256starstar.hpp"

#include "cxxopts.hpp"

#include <x86intrin.h>
#include <algorithm>
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

using weight_type = unsigned long;
using value_type = unsigned long;
static constexpr unsigned int retries = 100;

struct KnapsackInstance {
    struct Item {
        weight_type weight;
        value_type value;
    };
    std::vector<Item> items;
    weight_type capacity;
};

struct Node {
    Node() {}
    explicit Node(unsigned int i, weight_type cap, value_type v)
        : index{i}, free_capacity{cap}, value{v} {}
    unsigned int index;
    weight_type free_capacity;
    value_type value;
};

using payload_type = Node;

using PriorityQueue =
    typename multiqueue::MultiqueueFactory<value_type, payload_type,
                                           std::greater<>>::multiqueue_type;

struct Settings {
    std::filesystem::path instance_file;
    unsigned int num_threads = 4;
    std::uint64_t seed = 1;
    util::PriorityQueueParameters pq_params;
};

#ifdef COUNT_STATS
struct StatCounters {
    std::size_t pushed_nodes = 0;
    std::size_t ignored_nodes = 0;
    std::size_t extracted_nodes = 0;
    std::size_t failed_extracts = 0;
    std::size_t allocated_nodes = 0;
    std::size_t freed_nodes = 0;
    std::size_t processed_nodes = 0;
    std::size_t idle = 0;
};

static std::unique_ptr<StatCounters[]> stats;

inline void pushed_node(unsigned int id) { ++stats[id].pushed_nodes; }
inline void ignored_node(unsigned int id) { ++stats[id].ignored_nodes; }
inline void extracted_node(unsigned int id) { ++stats[id].extracted_nodes; }
inline void extract_failed(unsigned int id) { ++stats[id].failed_extracts; }
inline void allocated_node(unsigned int id) { ++stats[id].allocated_nodes; }
inline void freed_node(unsigned int id) { ++stats[id].freed_nodes; }
inline void processed_node(unsigned int id) { ++stats[id].processed_nodes; }
inline void started_idling(unsigned int id) { ++stats[id].idle; }
#else
inline void pushed_node(unsigned int) {}
inline void ignored_node(unsigned int) {}
inline void extracted_node(unsigned int) {}
inline void extract_failed(unsigned int) {}
inline void allocated_node(unsigned int) {}
inline void freed_node(unsigned int) {}
inline void processed_node(unsigned int) {}
inline void started_idling(unsigned int) {}
#endif

static Settings settings;
static KnapsackInstance instance;

value_type get_lower_bound(weight_type capacity, unsigned int index) noexcept {
    value_type value = 0;
    while (index < instance.items.size() &&
           instance.items[index].weight <= capacity) {
        capacity -= instance.items[index].weight;
        value += instance.items[index].value;
        ++index;
    }
    return value;
}

value_type get_upper_bound(weight_type capacity, unsigned int index) noexcept {
    value_type value = 0;
    while (index < instance.items.size() &&
           instance.items[index].weight <= capacity) {
        capacity -= instance.items[index].weight;
        value += instance.items[index].value;
        ++index;
    }
    if (index < instance.items.size() && capacity > 0) {
        value += static_cast<value_type>(
            instance.items[index].value * capacity /
                      instance.items[index].weight);
    }
    return value;
}

// Each thread has a state, which is either working (0), check_idle (1) or idle
// (2) or woken up by another thread (3)
struct alignas(2 * L1_CACHE_LINESIZE) IdleState {
    std::atomic_uint state = 0;
};
static std::unique_ptr<IdleState[]> idle_state;

alignas(2 * L1_CACHE_LINESIZE) std::atomic_size_t idle_counter = 0;

static inline bool idle(unsigned int id) {
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

// The best known solution value
alignas(2 * L1_CACHE_LINESIZE) std::atomic<value_type> best_value;

// Returns true if this node inserted new nodes
bool process_node(PriorityQueue::Handle& handle,
                  PriorityQueue::value_type const& value, unsigned int id) {
    Node const& node = value.second;
    if (node.index == instance.items.size()) {
        return false;
    }
    auto current_best_value = best_value.load(std::memory_order_relaxed);
    value_type upper_bound = value.first;
    if (upper_bound <= current_best_value) {
        // The upper bound of this node is worse than the currently best value
        ignored_node(id);
        return false;
    }
    processed_node(id);

    // Check if there is enough capacity for the next item
    if (instance.items[node.index].weight <= node.free_capacity) {
        value_type new_value = node.value + instance.items[node.index].value;
        while (new_value > current_best_value &&
               !best_value.compare_exchange_weak(current_best_value, new_value,
                                                 std::memory_order_relaxed)) {
            _mm_pause();
        }
        handle.push(
            {value.first,
             Node{node.index + 1,
                  node.free_capacity - instance.items[node.index].weight,
                  new_value}});
        pushed_node(id);
    }
    value_type upper_bound_without_next =
        node.value + get_upper_bound(node.free_capacity, node.index + 1);
    handle.push({upper_bound_without_next,
                 Node{node.index + 1, node.free_capacity, node.value}});
    pushed_node(id);
    return true;
}

void main_loop(typename PriorityQueue::Handle& handle, unsigned int id) {
    auto extract_item = [&handle, id](auto& retval) {
        for (std::size_t i = 0; i < retries * settings.num_threads; ++i) {
            if (handle.try_extract_top(retval)) {
                return true;
            }
            extract_failed(id);
        }
        return false;
    };
    auto all_empty = [&handle, id]() {
        for (unsigned int i = 0; i < *settings.pq_params.c; ++i) {
            if (!handle.is_empty(id * (*settings.pq_params.c) + i)) {
                return false;
            }
        }
        return true;
    };
    PriorityQueue::value_type retval;
    while (true) {
        if (!extract_item(retval)) {
            // no item found, initiate idling
            idle_state[id].state.store(1, std::memory_order_relaxed);
            idle_counter.fetch_add(1, std::memory_order_release);
            if (all_empty()) {
                if (idle(id)) {
                    return;
                }
            } else {
                idle_state[id].state.store(0, std::memory_order_relaxed);
                idle_counter.fetch_sub(1, std::memory_order_release);
            }
            continue;
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
            value_type upper_bound = get_upper_bound(instance.capacity, 0);
            value_type lower_bound = get_lower_bound(instance.capacity, 0);
            best_value.store(lower_bound, std::memory_order_relaxed);
            std::clog << "Solving knapsack instance..." << std::flush;
            if (lower_bound < upper_bound) {
                process_node(handle, {upper_bound, Node{0, instance.capacity, 0}},
                             ctx.get_id());
            }
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

static void read_problem() {
    std::ifstream file_stream{settings.instance_file};
    if (!file_stream) {
        throw std::runtime_error{"Could not open knapsack file"};
    }
    size_t n;
    file_stream >> n;
    if (!file_stream || file_stream.eof()) {
        throw std::runtime_error{"Error reading knapsack file"};
    }
    instance.items.clear();
    instance.items.reserve(n);
    file_stream >> instance.capacity;
    if (!file_stream || (n > 0 && file_stream.eof())) {
        throw std::runtime_error{"Error reading knapsack file"};
    }

    for (size_t i = 0; i < n; ++i) {
        if (!file_stream || file_stream.eof()) {
            throw std::runtime_error{"Error reading knapsack file"};
        }
        auto item = KnapsackInstance::Item{};
        file_stream >> item.value;
        if (!file_stream || file_stream.eof()) {
            throw std::runtime_error{"Error reading knapsack file"};
        }
        file_stream >> item.weight;
        if (!file_stream) {
            throw std::runtime_error{"Error reading knapsack file"};
        }
        instance.items.push_back(item);
    }
    std::sort(instance.items.begin(), instance.items.end(),
              [](auto const& lhs, auto const& rhs) {
                  return (static_cast<double>(lhs.value) /
                          static_cast<double>(lhs.weight)) >
                         (static_cast<double>(rhs.value) /
                          static_cast<double>(rhs.weight));
              });
}

int main(int argc, char* argv[]) {
    std::clog << "Build configuration:\n\n";
#ifndef NDEBUG
    std::clog << "DEBUG build\n";
#endif
    std::clog << "Payload: explicit\n";
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
        "Knapsack benchmark",
        "This executable measures and records the performance of relaxed "
        "priority queues in the knapsack problem");
    // clang-format off
    options.add_options()
      ("j,threads", "Specify the number of threads "
       "(default: 4)", cxxopts::value<unsigned int>(), "NUMBER")
      ("f,file", "The input graph", cxxopts::value<std::filesystem::path>(settings.instance_file)->default_value("knapsack.kp"), "PATH")
      ("s,seed", "Specify the initial seed"
       "(default: 0)", cxxopts::value<std::uint32_t>(), "NUMBER")
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
        if (result.count("threads") > 0) {
            settings.num_threads = result["threads"].as<unsigned int>();
        }
        if (result.count("seed") > 0) {
            settings.seed = result["seed"].as<std::uint32_t>();
        }
        if (result.count("factor") > 0) {
            settings.pq_params.c = result["factor"].as<std::size_t>();
        }
        if (result.count("stickiness") > 0) {
            settings.pq_params.stickiness =
                result["stickiness"].as<unsigned int>();
        }
    } catch (cxxopts::OptionParseException const& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    std::clog << "Settings: \n\t"
              << "Threads: " << settings.num_threads << "\n\t"
              << "Instance file: " << settings.instance_file.string() << "\n\t"
              << "Seed: " << settings.seed;
    std::clog << "\n\n";

    std::clog << "Reading problem..." << std::flush;
    try {
        read_problem();
    } catch (std::runtime_error const& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
    std::clog << "done\n\n";

    std::cout << "items: " << instance.items.size() << '\n';
    std::cout << "capacity: " << instance.capacity << '\n';

    idle_counter = 0;
    idle_state = std::make_unique<IdleState[]>(settings.num_threads);
#ifdef COUNT_STATS
    stats = std::make_unique<StatCounters[]>(settings.num_threads);
#endif

    xoroshiro256starstar rng;
    rng.seed(settings.seed);
    settings.pq_params.seed = rng();
    auto pq = util::create_pq<PriorityQueue>(1'000'000, settings.num_threads,
                                             settings.pq_params);

    std::clog << "Using priority queue: " << pq.description() << "\n\n";

    thread_coordination::ThreadCoordinator coordinator{settings.num_threads};
    coordinator.run_task<Task>(std::ref(pq));
    coordinator.join();
    auto duration = *coordinator.get_duration("main_loop");
    std::cout << "value: " << best_value << '\n';
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
            lhs.allocated_nodes += rhs.allocated_nodes;
            lhs.freed_nodes += rhs.freed_nodes;
            lhs.processed_nodes += rhs.processed_nodes;
            lhs.idle += rhs.idle;
            return lhs;
        });
    std::cout << "pushed nodes: " << total_stats.pushed_nodes << '\n';
    std::cout << "ignored nodes: " << total_stats.ignored_nodes << '\n';
    std::cout << "extracted nodes: " << total_stats.extracted_nodes << '\n';
    std::cout << "failed extracts: " << total_stats.failed_extracts << '\n';
    std::cout << "allocated nodes: " << total_stats.allocated_nodes << '\n';
    std::cout << "freed nodes: " << total_stats.freed_nodes << '\n';
    std::cout << "processed nodes: " << total_stats.processed_nodes << '\n';
    std::cout << "idle: " << total_stats.idle << '\n';
#endif
}
