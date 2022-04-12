#if !(defined PACKED_VALUE || defined HEAP_VALUE || defined EXPLICIT_VALUE)
#define PACKED_VALUE
#endif

#if !(defined LINEAR_MODE || defined BINARY_MODE || defined HINT_MODE)
#define HINT_MODE
#endif

#include "utils/priority_queue_factory.hpp"
#include "utils/thread_coordination.hpp"
#include "utils/threading.hpp"
#include "utils/xoroshiro256starstar.hpp"

#include "cxxopts.hpp"

#ifdef HEAP_VALUE
#include <tbb/scalable_allocator.h>
#endif
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

// Some pqs need higher values as sentinels
static constexpr value_type value_max =
    std::numeric_limits<value_type>::max() >> 2;

struct KnapsackInstance {
    struct Item {
        weight_type weight;
        value_type value;
    };
    std::vector<Item> items;
    weight_type capacity;
    std::vector<KnapsackInstance::Item> prefix_sum;
};

#ifdef HINT_MODE
struct Node {
    Node() {}
    explicit Node(unsigned int i, weight_type cap, value_type v, unsigned int h)
        : index{i}, free_capacity{cap}, value{v}, hint{h} {}
    unsigned int index;
    weight_type free_capacity;
    value_type value;
    unsigned int hint;
};
#else
struct Node {
    Node() {}
    explicit Node(unsigned int i, weight_type cap, value_type v)
        : index{i}, free_capacity{cap}, value{v} {}
    unsigned int index;
    weight_type free_capacity;
    value_type value;
};
#endif

#if defined PACKED_VALUE
using payload_type = unsigned long;
static_assert(sizeof(unsigned long) == sizeof(std::uint64_t),
              "64bit unsigned long required");

#ifdef HINT_MODE
static constexpr unsigned int bits_for_index = 12;
static constexpr unsigned int bits_for_free_capacity = 20;
static constexpr unsigned int bits_for_value = 20;
static constexpr unsigned int bits_for_hint = 12;
static_assert(bits_for_index + bits_for_hint + bits_for_free_capacity +
                      bits_for_value <=
                  64,
              "Can't use more than 64 bits to encode Node");
#else
static constexpr unsigned int bits_for_index = 16;
static constexpr unsigned int bits_for_free_capacity = 24;
static constexpr unsigned int bits_for_value = 24;
static_assert(bits_for_index + bits_for_free_capacity + bits_for_value <= 64,
              "Can't use more than 64 bits to encode Node");
#endif

static constexpr unsigned int index_pos = 0;
static constexpr unsigned int free_capacity_pos = index_pos + bits_for_index;
static constexpr unsigned int value_pos =
    free_capacity_pos + bits_for_free_capacity;
#ifdef HINT_MODE
static constexpr unsigned int hint_pos = value_pos + bits_for_value;
#endif

constexpr unsigned int get_index(payload_type p) noexcept {
    return static_cast<unsigned int>(
        (p >> index_pos) & p &
        ((static_cast<payload_type>(1) << bits_for_index) - 1));
}

constexpr weight_type get_free_capacity(payload_type p) noexcept {
    return static_cast<weight_type>(
        (p >> free_capacity_pos) &
        ((static_cast<payload_type>(1) << bits_for_free_capacity) - 1));
}

constexpr value_type get_value(payload_type p) noexcept {
    return static_cast<weight_type>(
        (p >> value_pos) &
        ((static_cast<payload_type>(1) << bits_for_value) - 1));
}

#ifdef HINT_MODE
constexpr unsigned int get_hint(payload_type p) noexcept {
    return static_cast<unsigned int>(
        (p >> (hint_pos)) &
        ((static_cast<payload_type>(1) << bits_for_hint) - 1));
}
#endif

constexpr payload_type to_payload(Node const& node) noexcept {
    return (static_cast<payload_type>(node.index) << index_pos) |
           (static_cast<payload_type>(node.free_capacity)
            << free_capacity_pos) |
           (static_cast<payload_type>(node.value) << value_pos)
#ifdef HINT_MODE
           | (static_cast<payload_type>(node.hint) << hint_pos)
#endif
        ;
}

#elif defined HEAP_VALUE

using payload_type = unsigned long;
static_assert(sizeof(payload_type) >= sizeof(std::uintptr_t),
              "Pointer must fit into value type of pq");

using Alloc = tbb::scalable_allocator<Node>;
using AllocTraits = std::allocator_traits<Alloc>;

Alloc global_alloc;

struct node_deleter {
    void operator()(Node* ptr) const {
        AllocTraits::destroy(global_alloc, ptr);
        AllocTraits::deallocate(global_alloc, ptr, 1);
    }
};

#elif defined EXPLICIT_VALUE
using payload_type = Node;
#endif

#ifdef PQ_MQ
using PriorityQueue =
    typename multiqueue::MultiqueueFactory<value_type, payload_type,
                                           std::greater<>>::multiqueue_type;
#else
using PriorityQueue =
    typename util::PriorityQueueFactory<value_type, payload_type>::type;
#endif

#if defined PQ_MQ || defined PQ_MF || defined PQ_TBB_Q
static constexpr bool reverse_priority = false;
#else
static constexpr bool reverse_priority = true;
#endif

using steady_type = std::chrono::steady_clock;

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
    std::size_t processed_nodes = 0;
    std::size_t idle = 0;
};

static std::unique_ptr<StatCounters[]> stats;

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

#ifdef LINEAR_MODE

value_type get_upper_bound(weight_type capacity, unsigned int index) noexcept {
    value_type value = 0;
    while (index < instance.items.size() &&
           instance.items[index].weight <= capacity) {
        capacity -= instance.items[index].weight;
        value += instance.items[index].value;
        ++index;
    }
    if (index < instance.items.size() && capacity > 0) {
        value +=
            static_cast<value_type>(instance.items[index].value * capacity /
                                    instance.items[index].weight);
    }
    return value;
}

#elif defined BINARY_MODE

value_type get_upper_bound(weight_type capacity, unsigned int index) noexcept {
    assert(index <= instance.items.size());
    assert(hint <= instance.prefix_sum.size());
    value_type value_offset = instance.prefix_sum[index].value;
    weight_type target_capacity = instance.prefix_sum[index].weight + capacity;
    auto it = std::upper_bound(
        instance.prefix_sum.begin(), instance.prefix_sum.end(), target_capacity,
        [](auto lhs, auto rhs) { return lhs < rhs.weight; });
    value_type result = (it - 1)->value - value_offset;
    if (it != instance.prefix_sum.end()) {
        double fraction =
            static_cast<double>(target_capacity - (it - 1)->weight) /
            static_cast<double>(it->weight - (it - 1)->weight);
        result += static_cast<value_type>(
            static_cast<double>(it->value - (it - 1)->value) * fraction);
    }
    return result;
}

#else

value_type get_upper_bound(weight_type capacity, unsigned int index,
                           unsigned int& hint) noexcept {
    assert(index <= instance.items.size());
    assert(hint <= instance.prefix_sum.size());
    value_type value_offset = instance.prefix_sum[index].value;
    weight_type target_capacity = instance.prefix_sum[index].weight + capacity;
    while (hint != instance.prefix_sum.size()) {
        if (instance.prefix_sum[hint].weight > target_capacity) {
            double fraction =
                static_cast<double>(target_capacity -
                                    instance.prefix_sum[hint - 1].weight) /
                static_cast<double>(instance.prefix_sum[hint].weight -
                                    instance.prefix_sum[hint - 1].weight);
            return (instance.prefix_sum[hint - 1].value - value_offset) +
                   static_cast<value_type>(
                       static_cast<double>(
                           instance.prefix_sum[hint].value -
                           instance.prefix_sum[hint - 1].value) *
                       fraction);
        }
        ++hint;
    }
    // All items fit
    return instance.prefix_sum[hint - 1].value - value_offset;
}
#endif

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
#ifdef PACKED_VALUE
    if (get_index(value.second) == instance.items.size()) {
#elif defined HEAP_VALUE
    auto node_ptr = std::unique_ptr<Node, node_deleter>(
        reinterpret_cast<Node*>(value.second));
    if (node_ptr->index == instance.items.size()) {
#elif defined EXPLICIT_VALUE
    if (value.second.index == instance.items.size()) {
#endif
        return false;
    }
    auto current_best_value = best_value.load(std::memory_order_relaxed);
    value_type upper_bound =
        reverse_priority ? value_max - value.first : value.first;
    if (upper_bound <= current_best_value) {
        // The upper bound of this node is worse than the currently best value
        ignored_node(id);
        return false;
    }
    processed_node(id);

#ifdef PACKED_VALUE
    Node node{get_index(value.second), get_free_capacity(value.second),
              get_value(value.second)
#ifdef HINT_MODE
                  ,
              get_hint(value.second)
#endif
    };
#elif defined HEAP_VALUE
    Node node = *node_ptr;
#elif defined EXPLICIT_VALUE
    Node const& node = value.second;
#endif

    bool inserted = false;
    // Check if there is enough capacity for the next item
#ifdef HINT_MODE
    assert(node.index < node.hint);
    assert(!!(instance.items[node.index].weight <= node.free_capacity) ==
           !!(node.index + 1 != node.hint));
    if (node.index + 1 != node.hint) {
#else
    if (instance.items[node.index].weight <= node.free_capacity) {
#endif
        value_type new_value = node.value + instance.prefix_sum[node.index + 1].value - instance.prefix_sum[node.index].value;
        value_type new_capacity = node.free_capacity - (instance.prefix_sum[node.index + 1].weight - instance.prefix_sum[node.index].weight);
        while (new_value > current_best_value &&
               !best_value.compare_exchange_weak(current_best_value, new_value,
                                                 std::memory_order_relaxed)) {
            _mm_pause();
        }
#ifdef PACKED_VALUE
        auto payload = to_payload(Node{
            node.index + 1,
            new_capacity, new_value
#ifdef HINT_MODE
            ,
            node.hint
#endif
        });
        handle.push({value.first, payload});
#elif defined HEAP_VALUE
        ++node_ptr->index;
        node_ptr->free_capacity = new_capacity;
        node_ptr->value = new_value;
        handle.push(value);
        node_ptr.release();
#elif defined EXPLICIT_VALUE
    handle.push(
        {value.first,
         Node{node.index + 1,
              new_capacity, new_value
#ifdef HINT_MODE
              ,
              node.hint
#endif
         }});
#endif
        pushed_node(id);
        inserted = true;
    }
#ifdef HINT_MODE
    unsigned int hint_without_next = node.hint;
    value_type upper_bound_without_next =
        node.value +
        get_upper_bound(node.free_capacity, node.index + 1, hint_without_next);
#else
    value_type upper_bound_without_next =
        node.value + get_upper_bound(node.free_capacity, node.index + 1);
#endif
    if (upper_bound_without_next <= current_best_value) {
        return inserted;
    }
    if (reverse_priority) {
        upper_bound_without_next = value_max - upper_bound_without_next;
    }
#ifdef PACKED_VALUE
    handle.push({upper_bound_without_next,
                 to_payload(Node{node.index + 1, node.free_capacity, node.value
#ifdef HINT_MODE
                                 ,
                                 hint_without_next
#endif
                 })});
#elif defined HEAP_VALUE
    if (node_ptr) {
        ++node_ptr->index;
#ifdef HINT_MODE
        node_ptr->hint = hint_without_next;
#endif
        handle.push({upper_bound_without_next, value.second});
        node_ptr.release();
    } else {
        auto ptr = AllocTraits::allocate(global_alloc, 1);
        AllocTraits::construct(global_alloc, ptr, node.index + 1,
                               node.free_capacity, node.value
#ifdef HINT_MODE
                               ,
                               hint_without_next
#endif
        );
        handle.push(
            {upper_bound_without_next, reinterpret_cast<payload_type>(ptr)});
    }
#elif defined EXPLICIT_VALUE
handle.push({upper_bound_without_next, Node{node.index + 1,
                               node.free_capacity, node.value
#ifdef HINT_MODE
                               ,
                               hint_without_next
#endif
                               }});
#endif
    pushed_node(id);
    return true;
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
#ifdef HINT_MODE
            unsigned int hint = 1;
            value_type upper_bound =
                get_upper_bound(instance.capacity, 0, hint);
#else
            value_type upper_bound = get_upper_bound(instance.capacity, 0);
#endif
            value_type lower_bound = get_lower_bound(instance.capacity, 0);
            best_value.store(lower_bound, std::memory_order_relaxed);
            std::clog << "Solving knapsack instance..." << std::flush;
            if (lower_bound < upper_bound) {
                if (reverse_priority) {
                    upper_bound = value_max - upper_bound;
                }
#ifdef PACKED_VALUE
                process_node(
                    handle,
                    {upper_bound, to_payload(Node{0, instance.capacity, 0
#ifdef HINT_MODE
                                                  ,
                                                  hint
#endif
                                  })},
                    ctx.get_id());
#elif defined HEAP_VALUE
                auto node_ptr = AllocTraits::allocate(global_alloc, 1);
                AllocTraits::construct(global_alloc, node_ptr, 0,
                                       instance.capacity, 0
#ifdef HINT_MODE
                                       ,
                                       hint
#endif
                );
                process_node(
                    handle,
                    {upper_bound, reinterpret_cast<payload_type>(node_ptr)},
                    ctx.get_id());
#elif defined EXPLICIT_VALUE
                process_node(handle,
                             {upper_bound, Node{0, instance.capacity, 0
#ifdef HINT_MODE
                                                ,
                                                hint
#endif
                                           }},
                             ctx.get_id());
#endif
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
    instance.prefix_sum.reserve(instance.items.size() + 1);
    instance.prefix_sum = {KnapsackInstance::Item{0, 0}};
    std::inclusive_scan(instance.items.begin(), instance.items.end(),
                        std::back_inserter(instance.prefix_sum),
                        [](auto const& lhs, auto const& rhs) {
                            return KnapsackInstance::Item{
                                lhs.weight + rhs.weight, lhs.value + rhs.value};
                        });
}

int main(int argc, char* argv[]) {
    std::clog << "Build configuration:\n\n";
#ifndef NDEBUG
    std::clog << "DEBUG build\n";
#endif
#ifdef PACKED_VALUE
    std::clog << "Payload: packed\n";
#elif defined HEAP_VALUE
    std::clog << "Payload: heap\n";
#else
    std::clog << "Payload: explicit\n";
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
#ifdef PACKED_VALUE
    if (instance.items.size() + 1 >= (1 << bits_for_index) ||
        instance.capacity >= (1 << bits_for_free_capacity) ||
        instance.prefix_sum.back().weight >= (1 << bits_for_free_capacity) ||
        instance.prefix_sum.back().value >= (1 << bits_for_value)) {
        std::cerr << "Instance cannot be represented in packed format\n";
        return 1;
    }
#endif

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
}
