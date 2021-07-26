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

using namespace std::chrono_literals;
using clock_type = std::chrono::steady_clock;

using weight_type = unsigned long;
using value_type = unsigned long;
using pointer_type = unsigned long;

static_assert(sizeof(pointer_type) >= sizeof(std::uintptr_t),
              "Pointer must fit into value type of pq");

// some pqs need sentinels
static constexpr value_type max_value =
    std::numeric_limits<value_type>::max() - 3;

using PriorityQueue =
    typename util::PriorityQueueFactory<value_type, pointer_type>::type;

static constexpr auto retries = 400;

struct Settings {
  std::filesystem::path instance_file;
  unsigned int num_threads = 4;
};

struct KnapsackInstance {
  struct Item {
    weight_type weight;
    value_type value;
  };
  std::vector<Item> items;
  weight_type capacity;
};

std::vector<std::vector<KnapsackInstance::Item>> prefix_sums;

struct Node {
  unsigned int index;
  weight_type weight;
  value_type value;
};

static value_type get_upper_bound(KnapsackInstance const& instance,
                                  Node const& node) {
  assert(node.index < instance.items.size());
  auto const free_capacity = instance.capacity - node.weight;
  auto last_to_include = std::lower_bound(
      prefix_sums[node.index].begin(), prefix_sums[node.index].end(),
      free_capacity, [&](auto const& prefix, auto const& capacity) {
        return prefix.weight < capacity;
      });
  if (last_to_include == prefix_sums[node.index].end()) {
    return node.value + prefix_sums[node.index].back().value;
  }
  if (last_to_include->weight == free_capacity) {
    return node.value + last_to_include->value;
  }
  auto const excess = last_to_include->weight - free_capacity;
  auto const& excess_item =
      instance.items[node.index +
                     (last_to_include - prefix_sums[node.index].begin())];
  assert(excess < excess_item.weight);
  return node.value + last_to_include->value -
         static_cast<value_type>(
             (static_cast<double>(excess) / excess_item.weight) *
             excess_item.value);
}

alignas(2 * L1_CACHE_LINESIZE) static std::atomic<value_type> best_value;

static std::atomic_size_t num_pushed_nodes = 0;
static std::atomic_size_t num_ignored_nodes = 0;
static std::atomic_size_t num_extracted_nodes = 0;
static std::atomic_size_t num_failed_extractions = 0;
static std::atomic_size_t num_allocated_nodes = 0;
static std::atomic_size_t num_freed_nodes = 0;
static std::atomic_size_t num_processed_nodes = 0;
static std::size_t* idle_count;

struct IdleState {
  alignas(2 * L1_CACHE_LINESIZE) std::atomic_uint state;
};

alignas(2 * L1_CACHE_LINESIZE) static std::atomic_size_t idle_counter = 0;
static IdleState* idle_state;

std::atomic_bool start_flag;

static inline bool idle(unsigned int id, unsigned int num_threads) {
  idle_state[id].state.store(2, std::memory_order_relaxed);
  idle_counter.fetch_add(1, std::memory_order_release);
  while (true) {
    if (idle_counter.load(std::memory_order_relaxed) == 2 * num_threads) {
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
                  KnapsackInstance const& instance) {
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
      Node* ptr = alloc_traits::allocate(alloc, 1);
      alloc_traits::construct(alloc, ptr, Node{0, 0, 0});
      ++num_local_allocated_nodes;
      // Some pqs need higher values as sentinels
      pq.push(handle, {max_value - get_upper_bound(instance, *ptr),
                       reinterpret_cast<pointer_type>(ptr)});
      ++num_local_pushed_nodes;
    }

    ctx.synchronize(stage++, [&ctx]() {
      std::clog << "Calculating knapsack solution...\n" << std::flush;
      ctx.notify_coordinator();
    });
    while (!start_flag.load(std::memory_order_relaxed)) {
      _mm_pause();
    }
    std::atomic_thread_fence(std::memory_order_acquire);
    std::pair<value_type, pointer_type> retval;
    while (true) {
      if (pq.extract_top(handle, retval)) {
      extracted:
        ++num_local_extracted_nodes;
        auto node_ptr = reinterpret_cast<Node*>(retval.second);
        auto current_node = *node_ptr;
        auto current_best_value = best_value.load(std::memory_order_relaxed);
        if (max_value - retval.first <= current_best_value) {
          ++num_local_ignored_nodes;
          alloc_traits::destroy(alloc, node_ptr);
          alloc_traits::deallocate(alloc, node_ptr, 1);
          ++num_local_freed_nodes;
          continue;
        }
        ++num_local_processed_nodes;
        if (current_node.index + 1 < instance.items.size()) {
          if (current_node.weight + instance.items[current_node.index].weight <=
              instance.capacity) {
            ++node_ptr->index;
            node_ptr->weight += instance.items[current_node.index].weight;
            node_ptr->value += instance.items[current_node.index].value;
            pq.push(handle, retval);
            node_ptr = nullptr;
            ++num_local_pushed_nodes;
          }
          ++current_node.index;
          auto const upper_bound_without_next =
              get_upper_bound(instance, current_node);
          if (upper_bound_without_next > current_best_value) {
            if (node_ptr == nullptr) {
              node_ptr = alloc_traits::allocate(alloc, 1);
              alloc_traits::construct(alloc, node_ptr, current_node);
              ++num_local_allocated_nodes;
            } else {
              ++node_ptr->index;
            }
            pq.push(handle, {max_value - upper_bound_without_next,
                             reinterpret_cast<pointer_type>(node_ptr)});
            node_ptr = nullptr;
            ++num_local_pushed_nodes;
          }
        } else {
          if (current_node.weight + instance.items[current_node.index].weight <=
              instance.capacity) {
            current_node.value += instance.items[current_node.index].value;
          }
          while (current_node.value > current_best_value &&
                 !best_value.compare_exchange_weak(current_best_value,
                                                   current_node.value)) {
          }
        }
        if (node_ptr == nullptr) {
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
          alloc_traits::destroy(alloc, node_ptr);
          alloc_traits::deallocate(alloc, node_ptr, 1);
          ++num_local_freed_nodes;
        }
      } else {
        ++num_local_failed_extractions;
        for (std::size_t i = 0; i < retries; ++i) {
          if (i < 20) {
            _mm_pause();
          } else {
            std::this_thread::yield();
          }
          if (pq.extract_top(handle, retval)) {
            goto extracted;
          }
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
  } catch (cxxopts::OptionParseException const& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }

#ifndef NDEBUG
  std::clog << "Using debug build!\n\n";
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
  std::clog << "\nItems\n"
            << std::setw(5) << "id" << std::setw(7) << "weight" << std::setw(7)
            << "value" << std::setw(6) << "ratio" << '\n';
  for (std::size_t i = 0; i < instance.items.size(); ++i) {
    std::clog << std::setw(5) << i << std::setw(8) << instance.items[i].weight
              << std::setw(8) << instance.items[i].value << std::setw(7)
              << std::setprecision(3)
              << static_cast<double>(instance.items[i].value) /
                     instance.items[i].weight
              << '\n';
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
  prefix_sums.resize(instance.items.size());
  for (std::size_t i = 0; i < instance.items.size(); ++i) {
    std::partial_sum(instance.items.begin() + i, instance.items.end(),
                     std::back_inserter(prefix_sums[i]),
                     [](auto const& lhs, auto const& rhs) {
                       return KnapsackInstance::Item{lhs.weight + rhs.weight,
                                                     lhs.value + rhs.value};
                     });
  }
  std::clog << "done\n";
#ifdef KNAPSACK_DEBUG_PRINT
  std::clog << "\nSorted items after ratio\n"
            << std::setw(5) << "id" << std::setw(7) << "weight" << std::setw(7)
            << "value" << std::setw(6) << "ratio" << '\n';
  for (std::size_t i = 0; i < instance.items.size(); ++i) {
    std::clog << std::setw(5) << i << std::setw(8) << instance.items[i].weight
              << std::setw(8) << instance.items[i].value << std::setw(7)
              << std::setprecision(3)
              << static_cast<double>(instance.items[i].value) /
                     instance.items[i].weight
              << '\n';
  }
  std::clog << "\nPrefix sum\n"
            << std::setw(5) << "id" << std::setw(7) << "weight" << std::setw(7)
            << "value" << '\n';
  for (std::size_t i = 0; i < prefix_sums[0].size(); ++i) {
    std::clog << std::setw(5) << i << std::setw(8) << prefix_sums[0][i].weight
              << std::setw(8) << prefix_sums[0][i].value << '\n';
  }
  std::clog << '\n';
#endif
  auto greedy_fit = std::lower_bound(
      prefix_sums[0].begin(), prefix_sums[0].end(), instance.capacity,
      [&](auto const& prefix, auto const& capacity) {
        return prefix.weight < capacity;
      });
  if (greedy_fit == prefix_sums[0].end()) {
    best_value = prefix_sums[0].back().value;
  } else {
    if (greedy_fit->value > instance.capacity) {
      best_value =
          (greedy_fit == prefix_sums[0].begin()) ? 0 : (greedy_fit - 1)->value;
    }
  }
  std::clog << "Greedy lower bound: " << best_value << '\n';
  idle_state = new IdleState[settings.num_threads]();

  idle_count = new std::size_t[settings.num_threads]();

  PriorityQueue pq{settings.num_threads};
  start_flag.store(false, std::memory_order_relaxed);
  std::atomic_thread_fence(std::memory_order_release);
  thread_coordination::ThreadCoordinator coordinator{settings.num_threads};
  coordinator.run<Task>(std::ref(pq), instance);
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
