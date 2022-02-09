#if !(defined PACKED_VALUE || defined HEAP_VALUE || defined EXPLICIT_VALUE)
#define PACKED_VALUE
#endif

#include "utils/operation_generator.hpp"
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

struct Node {
  Node() {}
  explicit Node(unsigned int i, unsigned int h, weight_type uc, value_type v)
      : index{i}, hint{h}, used_capacity{uc}, value{v} {}
  unsigned int index;
  unsigned int hint;
  weight_type used_capacity;
  value_type value;
};

#if defined PACKED_VALUE
using payload_type = unsigned long;
static_assert(sizeof(unsigned long) == sizeof(std::uint64_t),
              "64bit unsigned long required");

static constexpr unsigned int bits_for_index = 12;
static constexpr unsigned int bits_for_hint = 12;
static constexpr unsigned int bits_for_used_capacity = 20;
static constexpr unsigned int bits_for_value = 20;
static_assert(bits_for_index + bits_for_hint + bits_for_used_capacity +
                      bits_for_value <=
                  64,
              "Can't use more than 64 bits to encode Node");
constexpr unsigned int get_index(payload_type p) noexcept {
  return static_cast<unsigned int>(
      p & ((static_cast<payload_type>(1) << bits_for_index) - 1));
}

constexpr unsigned int get_hint(payload_type p) noexcept {
  return static_cast<unsigned int>(
      (p >> bits_for_index) &
      ((static_cast<payload_type>(1) << bits_for_hint) - 1));
}

constexpr weight_type get_used_capacity(payload_type p) noexcept {
  return static_cast<weight_type>(
      (p >> (bits_for_index + bits_for_hint)) &
      ((static_cast<payload_type>(1) << bits_for_used_capacity) - 1));
}

constexpr value_type get_value(payload_type p) noexcept {
  return static_cast<weight_type>(
      (p >> (bits_for_index + bits_for_hint + bits_for_used_capacity)) &
      ((static_cast<payload_type>(1) << bits_for_value) - 1));
}

constexpr payload_type to_payload(unsigned int index, unsigned int hint,
                                  weight_type used_capacity,
                                  value_type value) noexcept {
  return (static_cast<payload_type>(index) &
          ((static_cast<payload_type>(1) << bits_for_index) - 1)) |
         ((static_cast<payload_type>(hint) &
           ((static_cast<payload_type>(1) << bits_for_hint) - 1))
          << bits_for_index) |
         ((static_cast<payload_type>(used_capacity) &
           ((static_cast<payload_type>(1) << bits_for_used_capacity) - 1))
          << (bits_for_index + bits_for_hint)) |
         ((static_cast<payload_type>(value) &
           ((static_cast<payload_type>(1) << bits_for_value) - 1))
          << (bits_for_index + bits_for_hint + bits_for_used_capacity));
}

#elif defined HEAP_VALUE

using payload_type = unsigned long;
static_assert(sizeof(payload_type) >= sizeof(std::uintptr_t),
              "Pointer must fit into value type of pq");

using Alloc = tbb::scalable_allocator<Node>;
using AllocTraits = std::allocator_traits<Alloc>;

Alloc global_alloc;

#elif defined EXPLICIT_VALUE
using payload_type = Node;
#endif

#ifdef PQ_IS_MQ
using PriorityQueue = typename multiqueue::MultiqueueFactory<
    value_type, payload_type, std::greater<>>::template multiqueue_type<>;
#else
using PriorityQueue =
    typename util::PriorityQueueFactory<value_type, payload_type>::type;
#endif

using steady_type = std::chrono::steady_clock;

struct Settings {
  std::filesystem::path instance_file;
  unsigned int num_threads = 4;
  std::uint64_t seed = 1;
#if defined PQ_IS_MQ
  std::size_t c = PriorityQueue::param_type{}.c;
#ifndef PQ_MQ_RANDOM
  unsigned int stickiness = PriorityQueue::param_type{}.stickiness;
#endif
#endif
};

#ifdef COUNT_STATS
struct alignas(2 * L1_CACHE_LINESIZE) StatCounters {
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

value_type get_upper_bound(weight_type weight_unused, unsigned int index,
                           unsigned int& hint) noexcept {
  auto target_weight = instance.capacity + weight_unused;
  assert(hint <= instance.prefix_sum.size());
  while (hint != instance.prefix_sum.size() &&
         instance.prefix_sum[hint].weight <= target_weight) {
    ++hint;
  }
  value_type result =
      instance.prefix_sum[hint - 1].value - instance.prefix_sum[index].value;
  if (hint != instance.prefix_sum.size()) {
    weight_type free_capacity =
        target_weight - instance.prefix_sum[hint - 1].weight;
    result += (instance.items[hint - 1].value * free_capacity +
               instance.items[hint - 1].weight - 1) /
              instance.items[hint - 1].weight;
  }
  return result;
}

#ifdef EXACT_TERMINATION
// Each thread has a state, which is either working (0), check_idle (1) or idle
// (2) or woken up by another thread (3)
struct alignas(2 * L1_CACHE_LINESIZE) IdleState {
  std::atomic_uint state = 0;
};
static std::unique_ptr<IdleState[]> idle_state;

alignas(2 * L1_CACHE_LINESIZE) std::atomic_size_t idle_counter = 0;

static inline bool idle(unsigned int id) {
  idle_state[id].state.store(2, std::memory_order_relaxed);
  idle_counter.fetch_add(1, std::memory_order_release);
  started_idling(id);
  while (true) {
    if (idle_state[id].state.load(std::memory_order_relaxed) == 0) {
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
#endif

// The best known solution value
alignas(2 * L1_CACHE_LINESIZE) std::atomic<value_type> best_value;

// Returns true if this node was terminal
bool process_node(PriorityQueue::Handle& handle,
                  PriorityQueue::value_type const& value, unsigned int id) {
  auto current_best_value = best_value.load(std::memory_order_relaxed);
#ifdef PQ_IS_MQ
  value_type upper_bound = value.first;
#else
  value_type upper_bound = value_max - value.first;
#endif
  if (upper_bound <= current_best_value) {
    // The upper bound of this node is worse than the currently best value
    ignored_node(id);
#ifdef HEAP_VALUE
    AllocTraits::destroy(global_alloc, reinterpret_cast<Node*>(value.second));
    AllocTraits::deallocate(global_alloc, reinterpret_cast<Node*>(value.second),
                            1);
    freed_node(id);
#endif
    return true;
  }
  processed_node(id);

  bool terminal = true;
#ifdef PACKED_VALUE
  Node node{get_index(value.second), get_hint(value.second),
            get_used_capacity(value.second), get_value(value.second)};
#elif defined HEAP_VALUE
  Node* node_ptr = reinterpret_cast<Node*>(value.second);
  Node node = *node_ptr;
#elif defined EXPLICIT_VALUE
  Node const& node = value.second;
#endif
  while (node.value > current_best_value &&
         !best_value.compare_exchange_weak(current_best_value, node.value,
                                           std::memory_order_relaxed)) {
    _mm_pause();
  }
  // Check if there is enough capacity for the next item
  assert(node.index < node.hint);
  if (node.index + 1 != node.hint) {
#ifdef PACKED_VALUE
    auto payload =
        to_payload(node.index + 1, node.hint,
                   node.used_capacity + instance.items[node.index].weight,
                   node.value + instance.items[node.index].value);
    handle.push({value.first, payload});
#elif defined HEAP_VALUE
    node_ptr->index = node.index + 1;
    node_ptr->used_capacity =
        node.used_capacity + instance.items[node.index].weight;
    node_ptr->value = node.value + instance.items[node.index].value;
    handle.push(value);
#elif defined EXPLICIT_VALUE
    handle.push({value.first,
                 Node{node.index + 1, node.hint,
                      node.used_capacity + instance.items[node.index].weight,
                      node.value + instance.items[node.index].value}});
#endif
    terminal = false;
    pushed_node(id);
  }
  weight_type weight_unused =
      instance.prefix_sum[node.index + 1].weight - node.used_capacity;
  unsigned int hint_without_next = node.hint;
  value_type upper_bound_without_next =
      node.value +
      get_upper_bound(weight_unused, node.index + 1, hint_without_next);
  if (upper_bound_without_next <= current_best_value) {
#ifdef HEAP_VALUE
    if (terminal) {
      AllocTraits::destroy(global_alloc, node_ptr);
      AllocTraits::deallocate(global_alloc, node_ptr, 1);
      freed_node(id);
    }
#endif
    return terminal;
  }
  if (hint_without_next == instance.prefix_sum.size() ||
      instance.prefix_sum[hint_without_next - 1].weight ==
          instance.capacity + weight_unused) {
    // All remaining items fit or the upper bound matches
    while (upper_bound_without_next > current_best_value &&
           !best_value.compare_exchange_weak(current_best_value,
                                             upper_bound_without_next,
                                             std::memory_order_relaxed)) {
      _mm_pause();
    }
#ifdef HEAP_VALUE
    if (terminal) {
      AllocTraits::destroy(global_alloc, node_ptr);
      AllocTraits::deallocate(global_alloc, node_ptr, 1);
      freed_node(id);
    }
#endif
    return terminal;
  }
#ifndef PQ_IS_MQ
  upper_bound_without_next = value_max - upper_bound_without_next;
#endif
#ifdef PACKED_VALUE
  handle.push(
      {upper_bound_without_next, to_payload(node.index + 1, hint_without_next,
                                            node.used_capacity, node.value)});
#elif defined HEAP_VALUE
  if (terminal) {
    node_ptr->index = node.index + 1;
    node_ptr->hint = hint_without_next;
    handle.push({upper_bound_without_next, value.second});
  } else {
    node_ptr = AllocTraits::allocate(global_alloc, 1);
    AllocTraits::construct(global_alloc, node_ptr, node.index + 1,
                           hint_without_next, node.used_capacity, node.value);
    handle.push(
        {upper_bound_without_next, reinterpret_cast<payload_type>(node_ptr)});
    allocated_node(id);
  }
#elif defined EXPLICIT_VALUE
  handle.push({upper_bound_without_next, Node{node.index + 1, hint_without_next,
                                              node.used_capacity, node.value}});
#endif
  pushed_node(id);
  return false;
}

void main_loop(typename PriorityQueue::Handle& handle, unsigned int id) {
  PriorityQueue::value_type retval;
  while (true) {
    bool success = false;
    for (std::size_t i = 0; i < retries * settings.num_threads; ++i) {
      success = handle.try_extract_top(retval);
      if (success) {
        break;
      }
      extract_failed(id);
    }
    if (!success) {
#ifdef EXACT_TERMINATION
      idle_state[id].state.store(1, std::memory_order_relaxed);
      idle_counter.fetch_add(1, std::memory_order_release);
#ifdef PQ_IS_MQ
      for (unsigned int i = 0; i < settings.c; ++i) {
        if (handle.try_extract_from(id * settings.c + i, retval)) {
          success = true;
          break;
        }
      }
#else
      success = handle.try_extract_top(retval);
#endif
      if (!success) {
        if (idle(id)) {
          break;
        } else {
          continue;
        }
      }
      idle_state[id].state.store(0, std::memory_order_relaxed);
      idle_counter.fetch_sub(1, std::memory_order_release);
#else
      break;
#endif
    }
    extracted_node(id);
#ifdef EXACT_TERMINATION
    bool terminal = process_node(handle, retval, id);
    if (!terminal) {
      if (idle_counter.load(std::memory_order_acquire) != 0) {
        for (std::size_t i = 0; i < settings.num_threads; ++i) {
          if (i == id) {
            continue;
          }
          unsigned int thread_idle_state =
              idle_state[i].state.load(std::memory_order_relaxed);
          while (true) {
            while (thread_idle_state == 1) {
              std::this_thread::yield();
              thread_idle_state =
                  idle_state[i].state.load(std::memory_order_relaxed);
            }
            if (thread_idle_state != 2) {
              break;
            }
            if (idle_state[i].state.compare_exchange_strong(
                    thread_idle_state, 3, std::memory_order_relaxed)) {
              idle_counter.fetch_sub(2, std::memory_order_relaxed);
              idle_state[i].state.store(0, std::memory_order_release);
              break;
            }
          }
        }
      }
    }
#else
    process_node(handle, retval, id);
#endif
  }
}

struct Task {
  static void run(thread_coordination::Context ctx, PriorityQueue& pq) {
    auto handle = pq.get_handle();
    unsigned int hint;
    value_type upper_bound = get_upper_bound(0, 0, hint);
    value_type lower_bound = instance.prefix_sum[hint - 1].value;
    if (ctx.is_main()) {
      best_value.store(lower_bound, std::memory_order_relaxed);
      if (lower_bound == upper_bound) {
        std::clog << "Instance is trivial\n";
        return;
      }
      std::clog << "Solving knapsack instance..." << std::flush;
#ifndef PQ_IS_MQ
      upper_bound = value_max - upper_bound;
#endif
#ifdef PACKED_VALUE
      process_node(handle, {upper_bound, to_payload(0, hint, 0, 0)},
                   ctx.get_id());
#elif defined HEAP_VALUE
      auto node_ptr = AllocTraits::allocate(global_alloc, 1);
      AllocTraits::construct(global_alloc, node_ptr, 0, hint, 0, 0);
      allocated_node(ctx.get_id());
      process_node(handle,
                   {upper_bound, reinterpret_cast<payload_type>(node_ptr)},
                   ctx.get_id());
#elif defined EXPLICIT_VALUE
      process_node(handle, {upper_bound, Node{0, hint, 0, 0}}, ctx.get_id());
#endif
    } else if (lower_bound == upper_bound) {
      return;
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
  instance.prefix_sum = {KnapsackInstance::Item{0, 0}};
  instance.prefix_sum.reserve(instance.items.size() + 1);
  std::inclusive_scan(instance.items.begin(), instance.items.end(),
                      std::back_inserter(instance.prefix_sum),
                      [](auto const& lhs, auto const& rhs) {
                        return KnapsackInstance::Item{lhs.weight + rhs.weight,
                                                      lhs.value + rhs.value};
                      });
}

int main(int argc, char* argv[]) {
  std::clog << "Build configuration:\n\n";
#ifndef NDEBUG
  std::clog << "DEBUG build\n";
#endif
#ifdef EXACT_TERMINATION
  std::clog << "Termination: exact\n";
#else
  std::clog << "Termination: heuristic\n";
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
#if defined PQ_IS_MQ
      ("c,factor", "The number of queues per thread"
       "(default: 4)", cxxopts::value<std::size_t>(), "NUMBER")
#ifndef PQ_MQ_RANDOM
      ("k,stickiness", "The stickiness"
       "(default: 8)", cxxopts::value<unsigned int>(), "NUMBER")
#endif
#endif
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
#if defined PQ_IS_MQ
    if (result.count("factor") > 0) {
      settings.c = result["factor"].as<std::size_t>();
    }
#ifndef PQ_MQ_RANDOM
    if (result.count("stickiness") > 0) {
      settings.stickiness = result["stickiness"].as<unsigned int>();
    }
#endif
#endif
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

#ifdef EXACT_TERMINATION
  idle_counter = 0;
  idle_state = std::make_unique<IdleState[]>(settings.num_threads);
#endif
#ifdef COUNT_STATS
  stats = std::make_unique<StatCounters[]>(settings.num_threads);
#endif

  xoroshiro256starstar rng;
  rng.seed(settings.seed);
#if defined PQ_IS_MQ
  auto params = PriorityQueue::param_type{};
  params.seed = rng();
  params.c = settings.c;
#ifndef PQ_MQ_RANDOM
  params.stickiness = settings.stickiness;
#endif

  auto pq = PriorityQueue{1'000'000, settings.num_threads, params};
#else
  auto pq = PriorityQueue{1'000'000, settings.num_threads};
#endif

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
