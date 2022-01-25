#if !(defined PACKED_VALUE || defined POINTER_VALUE || defined EXPLICIT_VALUE)
#define PACKED_VALUE
#endif

#include "utils/operation_generator.hpp"
#include "utils/priority_queue_factory.hpp"
#include "utils/thread_coordination.hpp"
#include "utils/threading.hpp"
#include "utils/xoroshiro256starstar.hpp"

#include "cxxopts.hpp"

#ifdef POINTER_VALUE
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
static constexpr unsigned int retries = 400;

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

#if defined PACKED_VALUE
using payload = unsigned long;
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
constexpr unsigned int get_index(payload p) noexcept {
  return static_cast<unsigned int>(
      p & ((static_cast<payload>(1) << bits_for_index) - 1));
}

constexpr unsigned int get_hint(payload p) noexcept {
  return static_cast<unsigned int>(
      (p >> bits_for_index) & ((static_cast<payload>(1) << bits_for_hint) - 1));
}

constexpr weight_type get_used_capacity(payload p) noexcept {
  return static_cast<weight_type>(
      (p >> (bits_for_index + bits_for_hint)) &
      ((static_cast<payload>(1) << bits_for_used_capacity) - 1));
}

constexpr value_type get_value(payload p) noexcept {
  return static_cast<weight_type>(
      (p >> (bits_for_index + bits_for_hint + bits_for_used_capacity)) &
      ((static_cast<payload>(1) << bits_for_value) - 1));
}

constexpr payload to_payload(unsigned int index, unsigned int hint,
                             weight_type used_capacity,
                             value_type value) noexcept {
  return (static_cast<payload>(index) &
          ((static_cast<payload>(1) << bits_for_index) - 1)) |
         ((static_cast<payload>(hint) &
           ((static_cast<payload>(1) << bits_for_hint) - 1))
          << bits_for_index) |
         ((static_cast<payload>(used_capacity) &
           ((static_cast<payload>(1) << bits_for_used_capacity) - 1))
          << (bits_for_index + bits_for_hint)) |
         ((static_cast<payload>(value) &
           ((static_cast<payload>(1) << bits_for_value) - 1))
          << (bits_for_index + bits_for_hint + bits_for_used_capacity));
}

#elif defined POINTER_VALUE
struct Node {
  unsigned int index;
  unsigned int hint;
  weight_type used_capacity;
  value_type value;
};

using payload = unsigned long;
static_assert(sizeof(payload) >= sizeof(std::uintptr_t),
              "Pointer must fit into value type of pq");

constexpr unsigned int get_index(payload p) noexcept {
  return reinterpret_cast<Node*>(p)->index;
}

constexpr unsigned int get_hint(payload p) noexcept {
  return reinterpret_cast<Node*>(p)->hint;
}

constexpr weight_type get_used_capacity(payload p) noexcept {
  return reinterpret_cast<Node*>(p)->used_capacity;
}

constexpr value_type get_value(payload p) noexcept {
  return reinterpret_cast<Node*>(p)->value;
}

constexpr payload to_payload(Node const& node) noexcept {
  return reinterpret_cast<payload>(&node);
}

using Alloc = tbb::scalable_allocator<Node>;
using AllocTraits = std::allocator_traits<Alloc>;

Alloc alloc;

#elif defined EXPLICIT_VALUE
struct Node {
  unsigned int index;
  unsigned int hint;
  weight_type used_capacity;
  value_type value;
};

using payload = Node;

constexpr unsigned int get_index(payload const& p) noexcept { return p.index; }

constexpr unsigned int get_hint(payload const& p) noexcept { return p.hint; }

constexpr weight_type get_used_capacity(payload const& p) noexcept {
  return p.used_capacity;
}

constexpr value_type get_value(payload const& p) noexcept { return p.value; }

constexpr payload const& to_payload(Node const& node) noexcept { return node; }

#endif

#ifdef PQ_IS_MQ
using PriorityQueue = typename multiqueue::MultiqueueFactory<
    value_type, payload, std::greater<>>::template multiqueue_type<>;
#else
using PriorityQueue =
    typename util::PriorityQueueFactory<value_type, payload>::type;
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
  std::size_t failed_extractions = 0;
  std::size_t allocated_nodes = 0;
  std::size_t freed_nodes = 0;
  std::size_t processed_nodes = 0;
  std::size_t idle = 0;
};

static std::unique_ptr<StatCounters[]> stats;
#endif

static Settings settings;
static KnapsackInstance instance;

unsigned int get_upper_bound(weight_type used_capacity, unsigned int index,
                             unsigned int hint) noexcept {
  auto const weight_diff = instance.prefix_sum[index].weight - used_capacity;
  assert(hint < instance.prefix_sum.size());
  while (hint < instance.prefix_sum.size() &&
         instance.prefix_sum[hint].weight - weight_diff < instance.capacity) {
    ++hint;
  }
  return hint;
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
  while (true) {
    if (idle_state[id].state.load(std::memory_order_relaxed) == 0) {
      // Someone woke this thread up
      return false;
    }
    if (idle_counter.load(std::memory_order_relaxed) == 2 * num_threads) {
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
#ifdef COUNT_STATS
    ++stats[id].ignored_nodes;
    ++stats[id].freed_nodes;
#endif
#ifdef POINTER_VALUE
    AllocTraits::destroy(alloc, reinterpret_cast<Node*>(value.second));
    AllocTraits::deallocate(alloc, reinterpret_cast<Node*>(value.second), 1);
#endif
    return true;
  }
#ifdef COUNT_STATS
  ++stats[id].processed_nodes;
#endif

  bool terminal = true;
#ifdef POINTER_VALUE
  Node* node_ptr = reinterpret_cast<Node*>(value.second);
  Node node = *node_ptr;
#endif
  // Check if there is enough capacity for the next item
  if (get_index(value.second) != get_hint(value.second)) {
#ifdef PACKED_VALUE
    handle.push(
        {value.first,
         to_payload(get_index(value.second) + 1, get_hint(value.second),
                    get_used_capacity(value.second) +
                        instance.items[get_index(value.second)].weight,
                    get_value(value.second) +
                        instance.items[get_index(value.second)].value)});
#elif defined POINTER_VALUE
    ++node_ptr->index;
    node_ptr->used_capacity += instance.items[get_index(value.second)].weight;
    node_ptr->value += instance.items[get_index(value.second)].value;
    handle.push(value);
#elif defined EXPLICIT_VALUE
    handle.push({value.first,
                 {get_index(value.second) + 1, get_hint(value.second),
                  get_used_capacity(value.second) +
                      instance.items[get_index(value.second)].weight,
                  get_value(value.second) +
                      instance.items[get_index(value.second)].value}});
#endif
    terminal = false;
#ifdef COUNT_STATS
    ++stats[id].pushed_nodes;
#endif
  }
  auto const value_diff = instance.prefix_sum[get_index(value.second)].value -
                          get_value(value.second);
  auto const hint_without_next =
      get_upper_bound(get_used_capacity(value.second),
                      get_index(value.second) + 1, get_hint(value.second));
  auto const upper_bound_without_next =
      (hint_without_next != instance.prefix_sum.size()
           ? instance.prefix_sum[hint_without_next].value
           : instance.prefix_sum.back().value) -
      value_diff;
  if (hint_without_next == instance.prefix_sum.size() ||
      instance.prefix_sum[hint_without_next].weight -
              (instance.prefix_sum[get_index(value.second)].weight -
               get_used_capacity(value.second)) ==
          instance.capacity) {
    // All remaining items fit or the upper bound matches
    while (upper_bound_without_next > current_best_value &&
           !best_value.compare_exchange_weak(current_best_value,
                                             upper_bound_without_next,
                                             std::memory_order_relaxed)) {
      _mm_pause();
    }
    return terminal;
  }
  if (upper_bound_without_next > current_best_value) {
#ifndef PQ_IS_MQ
    upper_bound_without_next = value_max - upper_bound_without_next;
#endif
#ifdef PACKED_VALUE
    handle.push(
        {upper_bound_without_next,
         to_payload(get_index(value.second) + 1, hint_without_next,
                    get_used_capacity(value.second), get_value(value.second))});
#elif defined POINTER_VALUE
    if (terminal) {
      ++node_ptr->index;
      handle.push({upper_bound_without_next, value.second});
    } else {
      node_ptr = AllocTraits::allocate(alloc, 1);
      AllocTraits::construct(alloc, node_ptr, node.index + 1, hint_without_next,
                             node.used_capacity, node.value);
      handle.push(
          {upper_bound_without_next, reinterpret_cast<payload>(node_ptr)});
#ifdef COUNT_STATS
      ++stats[id].allocated_nodes;
#endif
    }
#elif defined EXPLICIT_VALUE
    handle.push({upper_bound_without_next,
                 {get_index(value.second) + 1, hint_without_next,
                  get_used_capacity(value.second), get_value(value.second)}});
#endif
    terminal = false;
#ifdef COUNT_STATS
    ++stats[id].pushed_nodes;
#endif
  }
#ifdef POINTER_VALUE
  if (terminal) {
    AllocTraits::destroy(alloc, reinterpret_cast<Node*>(value.second));
    AllocTraits::deallocate(alloc, reinterpret_cast<Node*>(value.second), 1);
#ifdef COUNT_STATS
    ++stats[id].freed_nodes;
#endif
  }
#endif
  return terminal;
}

void main_loop(typename PriorityQueue::Handle& handle, unsigned int id) {
  PriorityQueue::value_type retval;
  while (true) {
    bool success = false;
    for (std::size_t i = 0; i < retries; ++i) {
      success = handle.try_extract_top(retval);
      if (success) {
        break;
      }
#ifdef COUNT_STATS
      ++stats[id].failed_extractions;
#endif
      if (i > 100) {
        std::this_thread::yield();
      }
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
#ifdef COUNT_STATS
        ++stats[id].idle;
#endif
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
#ifdef COUNT_STATS
    ++stats[id].extracted_nodes;
#endif
#ifdef EXACT_TERMINATION
    bool is_leaf = process_node(handle, retval, id);
    if (!is_leaf) {
      if (idle_counter.load(std::memory_order_acquire) != 0) {
        for (std::size_t i = 0; i < ctx.get_num_threads(); ++i) {
          if (i == id) {
            continue;
          }
          unsigned int thread_idle_state =
              idle_state[i].load(std::memory_order_relaxed);
          while (true) {
            while (thread_idle_state == 1) {
              std::this_thread::yield();
              thread_idle_state = idle_state[i].load(std::memory_order_relaxed);
            }
            if (thread_idle_state != 2) {
              break;
            }
            assert(thread_idle_state == 0 || thread_idle_state == 3);
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

    auto hint = get_upper_bound(0, 0, 0);
    if (ctx.is_main()) {
      auto upper_bound = (hint == instance.prefix_sum.size())
                             ? instance.prefix_sum.back().value
                             : instance.prefix_sum[hint].value;
      if (hint == instance.prefix_sum.size() ||
          instance.prefix_sum[hint].weight == instance.capacity) {
        std::clog << "Instance is trivial\n";
        best_value = upper_bound;
        return;
      }
      std::clog << "Solving knapsack instance..." << std::flush;
#ifndef PQ_IS_MQ
      upper_bound = value_max - upper_bound;
#endif
#ifdef PACKED_VALUE
      process_node(handle, {upper_bound, to_payload(0, hint, 0, 0)},
                   ctx.get_id());
#elif defined POINTER_VALUE
      auto node_ptr = AllocTraits::allocate(alloc, 1);
      AllocTraits::construct(alloc, node_ptr, 0, hint, 0, 0);
#ifdef COUNT_STATS
      ++stats[ctx.get_id()].allocated_nodes;
#endif
      process_node(handle, {upper_bound, reinterpret_cast<payload>(node_ptr)},
                   ctx.get_id());
#elif defined EXPLICIT_VALUE
      process_node(handle, {upper_bound, {0, hint, 0, 0}}, ctx.get_id());
#endif
    } else {
      if (hint == instance.prefix_sum.size() ||
          instance.prefix_sum[hint].weight == instance.capacity) {
        return;
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
#elif defined POINTER_VALUE
  std::clog << "Payload: pointer\n";
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

#ifdef COUNT_STATS
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

  PriorityQueue pq{settings.num_threads, params};
#else
  PriorityQueue pq{settings.num_threads};
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
        lhs.failed_extractions += rhs.failed_extractions;
        lhs.allocated_nodes += rhs.allocated_nodes;
        lhs.freed_nodes += rhs.freed_nodes;
        lhs.processed_nodes += rhs.processed_nodes;
        lhs.idle += rhs.idle;
        return lhs;
      });
  std::cout << "pushed nodes: " << total_stats.pushed_nodes << '\n';
  std::cout << "ignored nodes: " << total_stats.ignored_nodes << '\n';
  std::cout << "extracted nodes: " << total_stats.extracted_nodes << '\n';
  std::cout << "failed extractions: " << total_stats.failed_extractions << '\n';
  std::cout << "allocated nodes: " << total_stats.allocated_nodes << '\n';
  std::cout << "freed nodes: " << total_stats.freed_nodes << '\n';
  std::cout << "processed nodes: " << total_stats.processed_nodes << '\n';
  std::cout << "idle: " << total_stats.idle << '\n';
#endif
}
