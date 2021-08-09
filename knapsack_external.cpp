#include "cxxopts.hpp"
#include "system_config.hpp"
#include "utils/priority_queue_factory.hpp"
#include "utils/thread_coordination.hpp"
#include "utils/threading.hpp"

#include <tbb/scalable_allocator.h>
#include <x86intrin.h>
#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <chrono>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <new>
#include <random>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

#ifdef KNAPSACK_DEBUG_PRINT
#include <mutex>
std::mutex write_mutex;
#endif

using namespace std::chrono_literals;
using clock_type = std::chrono::steady_clock;

using weight_type = unsigned long;
using value_type = unsigned long;
using pointer_type = unsigned long;

static_assert(sizeof(pointer_type) >= sizeof(std::uintptr_t),
              "Pointer must fit into value type of pq");

// some pqs need sentinels
constexpr value_type max_value = std::numeric_limits<value_type>::max() - 3;

using PriorityQueue =
    typename util::PriorityQueueFactory<value_type, pointer_type>::type;

// Number of retries before idle state is entered
constexpr auto retries = 400;

struct Settings {
  std::filesystem::path instance_file;
  unsigned int num_threads = 4;
  bool verbose = false;
};

struct KnapsackInstance {
  struct Item {
    weight_type weight;
    value_type value;
  };
  std::vector<Item> items;
  weight_type capacity;
  weight_type min_weight;
  std::vector<std::vector<KnapsackInstance::Item>> prefix_sums;
};

struct Node {
  std::size_t index;
  weight_type free_capacity;
  value_type value;
};

struct StatCounters {
  std::size_t pushed_nodes = 0;
  std::size_t ignored_nodes = 0;
  std::size_t extracted_nodes = 0;
  std::size_t failed_extractions = 0;
  std::size_t allocated_nodes = 0;
  std::size_t freed_nodes = 0;
  std::size_t processed_nodes = 0;
};

StatCounters& operator+=(StatCounters& lhs, StatCounters const& rhs) {
  lhs.pushed_nodes += rhs.pushed_nodes;
  lhs.ignored_nodes += rhs.ignored_nodes;
  lhs.extracted_nodes += rhs.extracted_nodes;
  lhs.failed_extractions += rhs.failed_extractions;
  lhs.allocated_nodes += rhs.allocated_nodes;
  lhs.freed_nodes += rhs.freed_nodes;
  lhs.processed_nodes += rhs.processed_nodes;
  return lhs;
}

StatCounters operator+(StatCounters lhs, StatCounters const& rhs) {
  lhs += rhs;
  return lhs;
}

// The best known solution value
alignas(2 * L1_CACHE_LINESIZE) std::atomic<value_type> best_value;

// Each thread has a state, which is either working (0), check_idle (1) or idle
// (2)
struct alignas(2 * L1_CACHE_LINESIZE) IdleState {
  std::atomic_uint state;
};
IdleState* idle_state;

alignas(2 * L1_CACHE_LINESIZE) std::atomic_size_t idle_counter = 0;

std::atomic_bool start_flag = false;

value_type get_upper_bound(KnapsackInstance const& instance, std::size_t index,
                           weight_type capacity) {
  assert(index < instance.items.size());
  auto it = std::lower_bound(instance.prefix_sums[index].begin(),
                             instance.prefix_sums[index].end(), capacity,
                             [&](auto const& sum, auto total_weight) {
                               return sum.weight < total_weight;
                             });
  if (it != instance.prefix_sums[index].end()) {
    // Imprecise but faster than fractional items
    return it->value;
  }
  return instance.prefix_sums[index].back();
}

// Returns true if new items were added by this node
bool process_node(PriorityQueue& pq, KnapsackInstance& instance,
                  std::pair<value_type, pointer_type> const& pq_item,
                  StatCounters& counters
#ifdef KNAPSACK_DEBUG_PRINT
                  ,
                  unsigned int id
#endif
) {
  ++counters.extracted_nodes;
  auto node_ptr = reinterpret_cast<Node*>(pq_item.second);
  auto current_best_value = best_value.load(std::memory_order_relaxed);
#ifdef KNAPSACK_DEBUG_PRINT
  if (verbose) {
    auto l = std::scoped_lock{write_mutex};
    std::clog << '[' << std::setw(3) << id
              << "] Current best value:  " << current_best_value << '\n';
    std::clog << '[' << std::setw(3) << id << "] Extracting node ("
              << node_ptr->index << ", " << node_ptr->weight << ", "
              << node_ptr->value << ") with upper bound "
              << max_value - pq_item.first << '\n';
  }
#endif
  if (max_value - pq_item.first <= current_best_value) {
    // The upper bound of this node is worse than the currently best value
    ++num_local_ignored_nodes;
#ifdef KNAPSACK_DEBUG_PRINT
    if (verbose) {
      auto l = std::scoped_lock{write_mutex};
      std::clog << '[' << std::setw(3) << id << "] Skipping node\n";
    }
#endif
    alloc_traits::destroy(alloc, node_ptr);
    alloc_traits::deallocate(alloc, node_ptr, 1);
    ++num_local_freed_nodes;
    return false;
  }
  ++counters.processed_nodes;
  if (node.index + 1 < instance.items) {
    // We only generate new nodes if the next item is not the last
    auto node = *node_ptr;
    if (instance.items[node.index].weight <= node.free_capacity) {
      ++node_ptr->index;
      node_ptr->free_capacity -= instance.items[node.index].weight;
      node_ptr->value += instance.items[node.index].value;
#ifdef KNAPSACK_DEBUG_PRINT
      if (verbose) {
        auto l = std::scoped_lock{write_mutex};
        std::clog << '[' << std::setw(3) << id
                  << "] Pushing node with next item (" << node_ptr->index
                  << ", " << node_ptr->weight << ", " << node_ptr->value
                  << ") and upper bound " << max_value - pq_item.first << '\n';
      }
#endif
      pq.push(handle,
              {pq_item.first, reinterpret_cast<pointer_type>(node_ptr)});
      node_ptr = nullptr;
      ++num_local_pushed_nodes;
    }
    ++node.index;
    auto const ub_no_next =
        node.value + get_upper_bound(instance, node.index, node.free_capacity);
#ifdef KNAPSACK_DEBUG_PRINT
    if (verbose) {
      auto l = std::scoped_lock{write_mutex};
      std::clog << '[' << std::setw(3) << id << "] Node without next item ("
                << node_ptr->index << ", " << node_ptr->weight << ", "
                << node_ptr->value << ") has upper bound " << ub_no_next
                << '\n';
    }
#endif
    if (ub_no_next > current_best_value) {
      if (node_ptr == nullptr) {
        node_ptr = alloc_traits::allocate(alloc, 1);
        alloc_traits::construct(alloc, node_ptr, node);
        ++num_local_allocated_nodes;
      } else {
        ++node_ptr->index;
      }
#ifdef KNAPSACK_DEBUG_PRINT
      if (verbose) {
        auto l = std::scoped_lock{write_mutex};
        std::clog << '[' << std::setw(3) << id
                  << "] Pushing node without next item" << '\n';
      }
#endif
      pq.push(handle, {max_value - ub_no_next,
                       reinterpret_cast<pointer_type>(node_ptr)});
      node_ptr = nullptr;
      ++num_local_pushed_nodes;
    }
    if (node_ptr == nullptr) {
      return true;
    }
    alloc_traits::destroy(alloc, node_ptr);
    alloc_traits::deallocate(alloc, node_ptr, 1);
    ++num_local_freed_nodes;
    return false;
  }

  // The next item is the last
  auto const final_value =
      node_ptr->value +
      (instance.items[node_ptr->index].weight <= node_ptr->free_capacity
           ? instance.items[node_ptr->index].value
           : 0);
#ifdef KNAPSACK_DEBUG_PRINT
  if (verbose) {
    auto l = std::scoped_lock{write_mutex};
    std::clog << '[' << std::setw(3) << id
              << "] Node has only one item left, final value: " << final_value
              << '\n';
    if (final_value > current_best_value) {
      std::clog << '[' << std::setw(3) << id << "] Updating best value\n";
    }
  }
#endif
  while (final_value > current_best_value &&
         !best_value.compare_exchange_weak(current_best_value, final_value)) {
  }
  alloc_traits::destroy(alloc, node_ptr);
  alloc_traits::deallocate(alloc, node_ptr, 1);
  ++num_local_freed_nodes;
  return false;
}

// Returns true if all work is done
bool idle(thread_coordination::Context const& ctx) {
  idle_state[ctx.get_id()].state.store(2, std::memory_order_relaxed);
  idle_counter.fetch_add(1, std::memory_order_release);
  while (true) {
    if (idle_counter.load(std::memory_order_relaxed) ==
        2 * ctx.get_num_threads()) {
      return true;
    }
    if (idle_state[id].state.load(std::memory_order_relaxed) == 0) {
      return false;
    }
    std::this_thread::yield();
  }
}

struct Task {
  static void run(thread_coordination::Context ctx, PriorityQueue& pq,
                  KnapsackInstance const& instance, bool verbose) {
#ifdef PQ_SPRAYLIST
    pq.init_thread(ctx.get_num_threads());
#endif
    unsigned int stage = 0;

    auto handle = pq.get_handle(ctx.get_id());

    tbb::scalable_allocator<Node> alloc;
    using alloc_traits = std::allocator_traits<decltype(alloc)>;

    std::size_t num_local_pushed_nodes = 0;
    std::size_t num_local_extracted_nodes = 0;
    std::size_t num_local_ignored_nodes = 0;
    std::size_t num_local_failed_extractions = 0;
    std::size_t num_local_allocated_nodes = 0;
    std::size_t num_local_freed_nodes = 0;
    std::size_t num_local_processed_nodes = 0;

    if (ctx.is_main()) {
      Node* node_ptr = alloc_traits::allocate(alloc, 1);
      alloc_traits::construct(alloc, node_ptr, Node{0, instance.capacity, 0});
      ++num_local_allocated_nodes;
      auto const ub = get_upper_bound(instance, *node_ptr);
#ifdef KNAPSACK_DEBUG_PRINT
      if (verbose) {
        auto l = std::scoped_lock{write_mutex};
        std::clog << '[' << std::setw(3) << ctx.get_id() << "] Pushing node ("
                  << node_ptr->index << ", " << node_ptr->weight << ", "
                  << node_ptr->value << ") with upper bound " << ub << '\n';
      }
#endif
      pq.push(handle, {max_value - ub, reinterpret_cast<pointer_type>(ptr)});
      ++num_local_pushed_nodes;
    }

    std::pair<value_type, pointer_type> retval;
    ctx.synchronize(stage++, [&ctx]() { ctx.notify_coordinator(); });
    while (!start_flag.load(std::memory_order_relaxed)) {
      _mm_pause();
    }
    std::atomic_thread_fence(std::memory_order_acquire);
    unsigned int try_count = 0;
    while (true) {
      if (pq.extract_top(handle, retval)) {
        bool new_work = process_node(pq, instance, retval, counters
#ifdef KNAPSACK_DEBUG_PRINT
                                     ,
                                     ctx.get_id()
#endif
        );
        if (new_work) {
          if (idle_counter.load(std::memory_order_acquire) > 0) {
            for (std::size_t i = 0; i < ctx.get_num_threads(); ++i) {
              if (i == ctx.get_id()) {
                continue;
              }
              unsigned int thread_state = 2;
              while (!idle_state[i].state.compare_exchange_weak(
                         thread_state, 3, std::memory_order_acq_rel,
                         std::memory_order_acquire) &&
                     thread_state != 0 && thread_state != 3) {
                thread_state = 2;
                std::this_thread::yield();
              }
              if (thread_state == 2) {
                idle_counter.fetch_sub(2, std::memory_order_relaxed);
                idle_state[i].state.store(0, std::memory_order_release);
              }
            }
          }
        } else {
#ifdef KNAPSACK_DEBUG_PRINT
          if (verbose) {
            auto l = std::scoped_lock{write_mutex};
            std::clog << '[' << std::setw(3) << ctx.get_id()
                      << "] Node has no children\n";
          }
#endif
        }
      } else {
        ++num_local_failed_extractions;
        ++try_count;
        if (try_count > 100) {
          std::this_thread::yield();
        }
        if (try_count < retries) {
          continue;
        }
        idle_state[ctx.get_id()].state.store(1, std::memory_order_relaxed);
        idle_counter.fetch_add(1, std::memory_order_release);
#ifdef PQ_IS_WRAPPER
        if (pq.extract_top(handle, retval)) {
#else
        if (pq.extract_from_partition(handle, retval)) {
#endif
          idle_counter.fetch_sub(1, std::memory_order_relaxed);
          idle_state[ctx.get_id()].state.store(0, std::memory_order_release);
          goto extracted;
        }
        ++idle_count[ctx.get_id()];
        if (idle(ctx.get_id(), ctx.get_num_threads())) {
          break;
        } else {
          continue;
        }
      }
    }
    num_allocated_nodes += num_local_allocated_nodes;
    num_extracted_nodes += num_local_extracted_nodes;
    num_failed_extractions += num_local_failed_extractions;
    num_freed_nodes += num_local_freed_nodes;
    num_ignored_nodes += num_local_ignored_nodes;
    num_processed_nodes += num_local_processed_nodes;
    num_pushed_nodes += num_local_pushed_nodes;
  }

  static threading::thread_config get_config(
      thread_coordination::Context const& ctx) {
    threading::thread_config config;
    config.cpu_set.reset();
    config.cpu_set.set(ctx.get_id());
    return config;
  }
};

static KnapsackInstance read_problem(Settings const& settings) {
  std::ifstream file_stream{settings.instance_file};
  if (!file_stream) {
    throw std::runtime_error{"Could not open knapsack file"};
  }
  size_t n;
  file_stream >> n;
  if (!file_stream || file_stream.eof()) {
    throw std::runtime_error{"Error reading knapsack file"};
  }

  KnapsackInstance instance;
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
  return instance;
}

int main(int argc, char* argv[]) {
  Settings settings{};

  cxxopts::Options options(
      "Knapsack benchmark",
      "This executable measures and records the performance of relaxed "
      "priority queues in the knapsack problem");
  // clang-format off
    options.add_options()
      ("j,threads", "Specify the number of threads "
       "(default: 4)", cxxopts::value<unsigned int>(), "NUMBER")
      ("f,file", "The input graph", cxxopts::value<std::filesystem::path>(settings.instance_file)->default_value("knapsack.kp"), "PATH")
      ("v,verbose", "Verbose output", cxxopts::value<bool>(), "BOOL")
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
    if (result.count("verbose") > 0) {
      settings.verbose = true;
    }
  } catch (cxxopts::OptionParseException const& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }

#ifndef NDEBUG
  std::clog << "Using debug build!\n\n";
#endif
#ifdef KNAPSACK_DEBUG_PRINT
  std::clog << "With debug print capability\n\n";
#endif
  std::clog << "Settings: \n\t"
            << "Threads: " << settings.num_threads << "\n\t"
            << "Instance file: " << settings.instance_file.string() << "\n\t";
  std::clog << "\n\n";

  std::clog << "Using priority queue: " << PriorityQueue::description() << '\n';
  KnapsackInstance instance;
  std::clog << "Reading problem..." << std::flush;
  try {
    instance = read_problem(settings);
  } catch (std::runtime_error const& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
  std::clog << "done\n";
  std::clog << "Num items: " << instance.items.size() << '\n';
  std::clog << "Capacity: " << instance.capacity << '\n';
#ifdef KNAPSACK_DEBUG_PRINT
  if (settings.verbose) {
    std::clog << "\nItems\n"
              << std::setw(5) << "id" << std::setw(8) << "weight"
              << std::setw(8) << "value" << std::setw(7) << "ratio" << '\n';
    for (std::size_t i = 0; i < instance.items.size(); ++i) {
      std::clog << std::setw(5) << i << std::setw(8) << instance.items[i].weight
                << std::setw(8) << instance.items[i].value << std::setw(7)
                << std::setprecision(3)
                << static_cast<double>(instance.items[i].value) /
                       instance.items[i].weight
                << '\n';
    }
  }
#endif
  std::clog << "Preprocessing...";
  std::sort(instance.items.begin(), instance.items.end(),
            [](auto const& lhs, auto const& rhs) {
              return (static_cast<double>(lhs.value) /
                      static_cast<double>(lhs.weight)) >
                     (static_cast<double>(rhs.value) /
                      static_cast<double>(rhs.weight));
            });
  instance.prefix_sums.resize(instance.items.size());
  for (std::size_t i = 0; i < instance.items.size(); ++i) {
    instance.prefix_sums[i].push_back({0, 0});
    std::partial_sum(instance.items.begin() + i, instance.items.end(),
                     std::back_inserter(instance.prefix_sums[i]),
                     [](auto const& lhs, auto const& rhs) {
                       return KnapsackInstance::Item{lhs.weight + rhs.weight,
                                                     lhs.value + rhs.value};
                     });
    auto end = std::upper_bound(
        instance.prefix_sums[i].begin(), instance.prefix_sums[i].end(),
        instance.capacity, [](auto capacity, auto const& prefix) {
          return capacity < prefix.weight;
        });
    instance.prefix_sums[i].erase(end, instance.prefix_sums[i].end());
  }
  std::clog << "done\n";
#ifdef KNAPSACK_DEBUG_PRINT
  if (settings.verbose) {
    std::clog << "\nSorted items after ratio\n"
              << std::setw(5) << "id" << std::setw(8) << "weight"
              << std::setw(8) << "value" << std::setw(7) << "ratio" << '\n';
    for (std::size_t i = 0; i < instance.items.size(); ++i) {
      std::clog << std::setw(5) << i << std::setw(8) << instance.items[i].weight
                << std::setw(8) << instance.items[i].value << std::setw(7)
                << std::setprecision(3)
                << static_cast<double>(instance.items[i].value) /
                       instance.items[i].weight
                << '\n';
    }
    std::clog << "\nPrefix sum\n"
              << std::setw(5) << "id" << std::setw(8) << "weight"
              << std::setw(8) << "value" << '\n';
    for (std::size_t i = 0; i < prefix_sums[0].size(); ++i) {
      std::clog << std::setw(5) << i << std::setw(8) << prefix_sums[0][i].weight
                << std::setw(8) << prefix_sums[0][i].value << '\n';
    }
    std::clog << '\n';
  }
#endif
  auto greedy_fit = std::upper_bound(prefix_sums[0].begin(),
                                     prefix_sums[0].end(), instance.capacity,
                                     [&](auto capacity, auto const& prefix) {
                                       return prefix.weight < capacity;
                                     });
  if (greedy_fit == prefix_sums[0].end()) {
    best_value = prefix_sums[0].back().value;
  } else {
    best_value =
        (greedy_fit == prefix_sums[0].begin()) ? 0 : (greedy_fit - 1)->value;
  }
  std::clog << "Greedy lower bound: " << best_value << '\n';
  idle_state = new IdleState[settings.num_threads]();

  idle_count = new std::size_t[settings.num_threads]();

  PriorityQueue pq{settings.num_threads};
  start_flag.store(false, std::memory_order_relaxed);
  std::atomic_thread_fence(std::memory_order_release);
  thread_coordination::ThreadCoordinator coordinator{settings.num_threads};
  std::clog << "Calculating knapsack solution...\n" << std::flush;
  coordinator.run<Task>(std::ref(pq), instance, settings.verbose);
  coordinator.wait_until_notified();
  start_flag.store(true, std::memory_order_release);
  auto start_tick = clock_type::now();
  __asm__ __volatile__("" ::: "memory");
  coordinator.join();
  __asm__ __volatile__("" ::: "memory");
  auto end_tick = clock_type::now();
  std::clog << "Solution value: " << best_value << '\n';
  std::clog << "Time (ms): "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end_tick -
                                                                     start_tick)
                   .count()
            << '\n';
  std::clog << "Pushed nodes: " << num_pushed_nodes << '\n';
  std::clog << "Extracted nodes: " << num_extracted_nodes << '\n';
  std::clog << "Ignored nodes: " << num_ignored_nodes << '\n';
  std::clog << "Failed extractions: " << num_failed_extractions << '\n';
  std::clog << "Allocated nodes: " << num_allocated_nodes << '\n';
  std::clog << "Freed nodes: " << num_freed_nodes << '\n';
  std::clog << "Processed nodes: " << num_processed_nodes << '\n';
  std::clog << "Idle count per thread:\n";
  for (unsigned int t = 0; t < settings.num_threads; ++t) {
    std::clog << std::setw(12) << idle_count[t];
  }
  std::clog << '\n';
  std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end_tick -
                                                                     start_tick)
                   .count()
            << ' ' << num_processed_nodes << '\n';

  delete[] idle_state;
  delete[] idle_count;
  return 0;
}
