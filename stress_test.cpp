#include "external/cxxopts.hpp"

#include "system_config.hpp"
#include "utils/operation_generator.hpp"
#include "utils/priority_queue_factory.hpp"
#include "utils/thread_coordination.hpp"
#include "utils/threading.hpp"

#include <time.h>
#include <x86intrin.h>
#include <array>
#include <atomic>
#include <cassert>
#include <chrono>
#include <fstream>
#include <functional>
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

#if !defined THROUGHPUT && !defined QUALITY
#error Need to define either THROUGHPUT or QUALITY
#endif

#if defined THROUGHPUT
#undef QUALITY
#endif
#if defined QUALITY
#undef THROUGHPUT
#endif

#ifndef L1_CACHE_LINESIZE
#error Need to define L1_CACHE_LINESIZE
#endif

using key_type = unsigned long;
using value_type = unsigned long;

using PriorityQueue =
    typename util::PriorityQueueFactory<key_type, value_type>::type;

using namespace std::chrono_literals;

#ifdef QUALITY

using tick_type = std::uint64_t;

/* using clock_type = std::chrono::steady_clock; */
/* using time_point_type = clock_type::time_point; */
/* static inline time_point_type get_time_point() noexcept { */
/*   return clock_type::now(); */
/* } */

/* static inline tick_type get_tick() noexcept { */
/*   return
 * static_cast<tick_type>(clock_type::now().time_since_epoch().count()); */
/* } */

static inline tick_type get_tick() noexcept {
#ifdef USE_TSC
  // Not synchronized among sockets
  return __rdtsc();
#else
  timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);
  return static_cast<tick_type>(ts.tv_sec * 1000000000 + ts.tv_nsec);
#endif
}

#endif

struct Settings {
  std::size_t prefill_size = 1'000'000;
  std::size_t num_operations = 10'000'000;
#ifdef QUALITY
  std::chrono::microseconds sleep_between_operations =
      std::chrono::microseconds::zero();
#else
  unsigned int iterations = 5;
#endif
  unsigned int num_threads = 4;
  std::uint32_t seed = 1;
#if defined PQ_MQ_RANDOM || defined PQ_MQ_STICKY
  std::size_t c = 4;
#endif
#ifdef PQ_MQ_STICKY
  unsigned int stickiness = 8;
#endif
  OperationGenerator<key_type> insert_config{
      InsertPolicy::Uniform,
      KeyDistribution::Uniform,
      std::numeric_limits<value_type>::min(),
      static_cast<value_type>(std::numeric_limits<std::uint32_t>::max() -
                              1),  // Some pqs use sentinels
      1,
      100,
  };
};

struct OperationCount {
  std::size_t num_prefill_insertions = 0;
  std::size_t num_insertions = 0;
  std::size_t num_deletions = 0;
  std::size_t num_failed_deletions = 0;
};

OperationCount& operator+=(OperationCount& lhs,
                           OperationCount const& rhs) noexcept {
  lhs.num_prefill_insertions += rhs.num_prefill_insertions;
  lhs.num_insertions += rhs.num_insertions;
  lhs.num_deletions += rhs.num_deletions;
  lhs.num_failed_deletions += rhs.num_failed_deletions;
  return lhs;
}

OperationCount operator+(OperationCount lhs,
                         OperationCount const& rhs) noexcept {
  return lhs += rhs;
}

#ifdef QUALITY

static constexpr unsigned int bits_for_thread_id = 8;
static constexpr value_type value_mask =
    (static_cast<value_type>(1)
     << (std::numeric_limits<value_type>::digits - bits_for_thread_id)) -
    1;

static constexpr value_type to_value(unsigned int thread_id,
                                     value_type elem_id) noexcept {
  return (static_cast<value_type>(thread_id)
          << (std::numeric_limits<value_type>::digits - bits_for_thread_id)) |
         (elem_id & value_mask);
}

static constexpr unsigned int get_thread_id(value_type value) noexcept {
  return static_cast<unsigned int>(
      value >> (std::numeric_limits<value_type>::digits - bits_for_thread_id));
}

static constexpr value_type get_elem_id(value_type value) noexcept {
  return value & value_mask;
}

struct InsertionLogEntry {
  tick_type tick;
  key_type key;
};

struct DeletionLogEntry {
  tick_type tick;
  value_type value;
};

#endif

struct alignas(L1_CACHE_LINESIZE) ThreadData {
  OperationCount op_count;
#ifdef QUALITY
  std::vector<InsertionLogEntry> ins_log;
  std::vector<DeletionLogEntry> del_log;
#endif
  std::uint32_t seed;
  InsertingStrategy<key_type> inserter;
};

static Settings settings;
alignas(L1_CACHE_LINESIZE) static key_type* prefill_keys;
alignas(L1_CACHE_LINESIZE) static key_type* keys;
alignas(L1_CACHE_LINESIZE) static std::atomic_size_t key_index;
alignas(L1_CACHE_LINESIZE) static ThreadData* thread_data;

template <typename F, typename... Args>
void execute_blockwise(std::size_t n, F f, Args const&... args) {
  static constexpr std::size_t block_size = 4096;
  std::size_t current_begin, current_end;
  current_begin = key_index.load(std::memory_order_relaxed);
  while (current_begin < n) {
    current_end = std::min(current_begin + block_size, n);
    if (!key_index.compare_exchange_weak(current_begin, current_end,
                                         std::memory_order_relaxed)) {
      continue;
    }
    f(current_begin, current_end, args...);
  }
}

void generate_prefill_keys(unsigned int id) {
  execute_blockwise(
      settings.prefill_size, [id](std::size_t begin, std::size_t end) {
        std::generate(prefill_keys + begin, prefill_keys + end,
                      [id]() { return thread_data[id].inserter.get_key(); });
      });
}

void generate_keys(unsigned int id) {
  execute_blockwise(settings.num_operations,
                    [id](std::size_t begin, std::size_t end) {
                      std::generate(keys + begin, keys + end, [id]() {
                        return thread_data[id].inserter.insert()
                                   ? thread_data[id].inserter.get_key()
                                   : std::numeric_limits<std::uint32_t>::max();
                      });
                    });
}

#ifdef QUALITY
void prefill(unsigned int id, PriorityQueue::Handle& handle,
             PriorityQueue& pq) {
  execute_blockwise(settings.prefill_size, [id, &handle, &pq](std::size_t begin,
                                                              std::size_t end) {
    for (key_type* k = prefill_keys + begin; k != prefill_keys + end; ++k) {
      pq.push(
          handle,
          {*k,
           to_value(id, static_cast<value_type>(
                            thread_data[id].op_count.num_prefill_insertions))});
      thread_data[id].ins_log.push_back(InsertionLogEntry{0, *k});
      ++thread_data[id].op_count.num_prefill_insertions;
    }
  });
}

void work(unsigned int id, PriorityQueue::Handle& handle, PriorityQueue& pq) {
  std::seed_seq seq{thread_data[id].seed + 2};
  auto gen = std::mt19937(seq);
  auto dist = std::uniform_int_distribution<long>(
      0, settings.sleep_between_operations.count());
  execute_blockwise(settings.num_operations, [&, id](std::size_t begin,
                                                     std::size_t end) {
    for (key_type* k = keys + begin; k != keys + end; ++k) {
      if (*k == std::numeric_limits<std::uint32_t>::max()) {
        PriorityQueue::value_type retval;
        while (!pq.try_delete_min(handle, retval)) {
          // Not expected due to prefill
          ++thread_data[id].op_count.num_failed_deletions;
        }
        auto tick = get_tick();
#ifdef USE_TSC
        // https://stackoverflow.com/questions/54690703/solution-to-rdtsc-out-of-order-execution
        _mm_lfence();
#endif
        thread_data[id].del_log.push_back(DeletionLogEntry{tick, retval.data});
        ++thread_data[id].op_count.num_deletions;
      } else {
        value_type value =
            to_value(id, thread_data[id].op_count.num_prefill_insertions +
                             thread_data[id].op_count.num_insertions);
        pq.push(handle, {*k, value});
        auto tick = get_tick();
#ifdef USE_TSC
        _mm_lfence();
#endif
        thread_data[id].ins_log.push_back(InsertionLogEntry{tick, *k});
        ++thread_data[id].op_count.num_insertions;
      }
      if (settings.sleep_between_operations !=
          std::chrono::microseconds::zero()) {
        std::this_thread::sleep_for(std::chrono::microseconds{dist(gen)});
      }
    }
  });
}

#else

void prefill(unsigned int id, PriorityQueue::Handle& handle,
             PriorityQueue& pq) {
  execute_blockwise(
      settings.prefill_size, [&](std::size_t begin, std::size_t end) {
        for (key_type* k = prefill_keys + begin; k != prefill_keys + end; ++k) {
          pq.push(handle, {*k, *k});
          ++thread_data[id].op_count.num_insertions;
        }
      });
}

void work(unsigned int id, PriorityQueue::Handle& handle, PriorityQueue& pq) {
  PriorityQueue::value_type retval;
  // Let retval escape
  asm volatile("" ::"g"(&retval));
  execute_blockwise(settings.num_operations,
                    [&, id](std::size_t begin, std::size_t end) {
                      for (key_type* k = keys + begin; k != keys + end; ++k) {
                        if (*k == std::numeric_limits<std::uint32_t>::max()) {
                          while (!pq.try_delete_min(handle, retval)) {
                            // Not expected due to prefill
                            ++thread_data[id].op_count.num_failed_deletions;
                          }
                          // "Use" memory to force write to retval
                          asm volatile("" ::: "memory");
                          ++thread_data[id].op_count.num_deletions;
                        } else {
                          pq.push(handle, {*k, *k});
                          ++thread_data[id].op_count.num_insertions;
                        }
                      }
                    });
}
#endif

struct Task {
  static void run(thread_coordination::Context ctx, PriorityQueue& pq,
                  std::chrono::milliseconds& duration) {
#ifdef PQ_SPRAYLIST
    pq.init_thread(ctx.get_num_threads());
#endif

    unsigned int stage = 0;

    auto handle = pq.get_handle();

#ifdef PQ_IS_MQ
    handle.data.seed(thread_data[ctx.get_id()].seed);
#endif

    if (ctx.is_main()) {
      std::clog << "Generating operations..." << std::flush;
      key_index = 0;
    }
    ctx.execute_synchronized(stage++, generate_prefill_keys, ctx.get_id());
    ctx.synchronize(stage++, []() { key_index = 0; });
    ctx.execute_synchronized(stage++, generate_keys, ctx.get_id());
    ctx.synchronize(stage++, []() {
      std::clog << "done\nPrefilling..." << std::flush;
      key_index = 0;
    });
    ctx.execute_synchronized(stage++, prefill, ctx.get_id(), handle, pq);

    ctx.synchronize(stage++, [&ctx]() {
      std::clog << "done\nStarting the stress test..." << std::flush;
      key_index = 0;
    });
    ctx.synchronize(stage++);
    if (ctx.is_main()) {
      auto start_time = std::chrono::steady_clock::now();
      work(ctx.get_id(), handle, pq);
      ctx.synchronize(stage++);
      duration = std::chrono::duration_cast<std::chrono::milliseconds>(
          std::chrono::steady_clock::now() - start_time);
      std::clog << "done\n";
    } else {
      work(ctx.get_id(), handle, pq);
      ctx.synchronize(stage++);
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

int main(int argc, char* argv[]) {
#ifndef NDEBUG
  std::clog << "DEBUG build!\n\n";
#endif

  std::clog << "Build configuration:\n\t";
#if defined THROUGHPUT
  std::clog << "Mode: Throughput\n\t";
#elif defined QUALITY
  std::clog << "Mode: Quality\n\t";
#endif
  std::clog << "L1 cache linesize (byte): " << L1_CACHE_LINESIZE << "\n\t";
  std::clog << '\n';

  std::clog << "Command line: ";
  std::copy(argv, argv + argc,
            std::ostream_iterator<char const*>(std::clog, " "));
  std::clog << "\n\n";

  cxxopts::Options options(argv[0]);
  // clang-format off
    options.add_options()
      ("p,prefill", "Specify the number of elements to prefill the queue with "
       "(default: 1'000'000)", cxxopts::value<size_t>(), "NUMBER")
      ("n,ops", "Specify the number of operations "
       "(default: 10'000'000)", cxxopts::value<std::size_t>(), "NUMBER")
      ("i,insert", "Specify the insert policy as one of \"uniform\", \"split\", \"producer\", \"alternating\" "
       "(default: uniform)", cxxopts::value<std::string>(), "ARG")
      ("j,threads", "Specify the number of threads "
       "(default: 4)", cxxopts::value<unsigned int>(), "NUMBER")
  #ifdef QUALITY
      ("w,sleep", "Specify the sleep time between operations in microseconds"
       "(default: 0)", cxxopts::value<unsigned int>(), "NUMBER")
  #endif
      ("d,distribution", "Specify the key distribution as one of \"uniform\", \"dijkstra\", \"ascending\", \"descending\", \"threadid\" "
       "(default: uniform)", cxxopts::value<std::string>(), "ARG")
      ("m,max", "Specify the max key "
       "(default: MAX)", cxxopts::value<key_type>(), "NUMBER")
      ("l,min", "Specify the min key "
       "(default: 0)", cxxopts::value<key_type>(), "NUMBER")
      ("s,seed", "Specify the initial seed"
       "(default: 0)", cxxopts::value<std::uint32_t>(), "NUMBER")
  #if defined PQ_MQ_RANDOM || defined PQ_MQ_STICKY
      ("c,factor", "The number of queues per thread"
       "(default: 4)", cxxopts::value<std::size_t>(), "NUMBER")
  #endif
  #ifdef PQ_MQ_STICKY
      ("k,stickiness", "The stickiness"
       "(default: 8)", cxxopts::value<unsigned int>(), "NUMBER")
  #endif
      ("h,help", "Print this help");
  // clang-format on

  try {
    auto result = options.parse(argc, argv);
    if (result.count("help")) {
      std::cerr << options.help() << std::endl;
      return 0;
    }
    if (result.count("prefill") > 0) {
      settings.prefill_size = result["prefill"].as<size_t>();
    }
    if (result.count("ops") > 0) {
      settings.num_operations = result["ops"].as<std::size_t>();
    }
    if (result.count("insert") > 0) {
      std::string policy = result["insert"].as<std::string>();
      if (policy == "uniform") {
        settings.insert_config.insert_policy = InsertPolicy::Uniform;
      } else if (policy == "alternating") {
        settings.insert_config.insert_policy = InsertPolicy::Alternating;
      } else {
        std::cerr << "Unknown insert policy \"" << policy << "\"\n";
        return 1;
      }
    }
    if (result.count("distribution") > 0) {
      std::string dist = result["distribution"].as<std::string>();
      if (dist == "uniform") {
        settings.insert_config.key_distribution = KeyDistribution::Uniform;
      } else if (dist == "ascending") {
        settings.insert_config.key_distribution = KeyDistribution::Ascending;
      } else if (dist == "descending") {
        settings.insert_config.key_distribution = KeyDistribution::Descending;
      } else if (dist == "dijkstra") {
        settings.insert_config.key_distribution = KeyDistribution::Dijkstra;
      } else {
        std::cerr << "Unknown key distribution \"" << dist << "\"\n";
        return 1;
      }
    }
    if (result.count("threads") > 0) {
      settings.num_threads = result["threads"].as<unsigned int>();
    }
#ifdef QUALITY
    if (result.count("sleep") > 0) {
      settings.sleep_between_operations =
          std::chrono::microseconds{result["sleep"].as<unsigned int>()};
    }
#endif
    if (result.count("max") > 0) {
      settings.insert_config.max_key = result["max"].as<key_type>();
    }
    if (result.count("min") > 0) {
      settings.insert_config.min_key = result["min"].as<key_type>();
    }
    if (result.count("seed") > 0) {
      settings.seed = result["seed"].as<std::uint32_t>();
    }
#if defined PQ_MQ_RANDOM || defined PQ_MQ_STICKY
    if (result.count("factor") > 0) {
      settings.c = result["factor"].as<std::size_t>();
    }
#endif
#ifdef PQ_MQ_STICKY
    if (result.count("stickiness") > 0) {
      settings.stickiness = result["stickiness"].as<unsigned int>();
    }
#endif
  } catch (cxxopts::OptionParseException const& e) {
    std::cerr << "Error parsing arguments: " << e.what() << '\n';
    std::cerr << options.help() << std::endl;
    return 1;
  }

  std::clog << "Settings: \n\t"
            << "Prefill: " << settings.prefill_size << "\n\t"
            << "Operations: " << settings.num_operations << "\n\t"
            << "Threads: " << settings.num_threads << "\n\t"
            << "Insert policy: "
            << get_insert_policy_name(settings.insert_config.insert_policy)
            << "\n\t"
            << "Key distribution: "
            << get_key_distribution_name(
                   settings.insert_config.key_distribution)
            << "\n\t"
            << "Min key: " << settings.insert_config.min_key << "\n\t"
            << "Max key: " << settings.insert_config.max_key << "\n\t"
#ifdef QUALITY
            << "Sleep between operations: "
            << settings.sleep_between_operations.count() << " us\n\t"
#endif
            << "Dijkstra min increase: "
            << settings.insert_config.dijkstra_min_increase << "\n\t"
            << "Dijkstra max increase: "
            << settings.insert_config.dijkstra_max_increase << "\n\t"
            << "Seed: " << settings.seed;
  std::clog << "\n\n";

#ifdef QUALITY
  if (settings.num_threads > (1 << bits_for_thread_id) - 1) {
    std::cerr << "Too many threads, increase the number of thread bits!"
              << std::endl;
    return 1;
  }
#endif

#if defined PQ_MQ_RANDOM
  PriorityQueue pq{settings.num_threads, settings.c};
#elif defined PQ_MQ_STICKY
  PriorityQueue pq{settings.num_threads, settings.c, settings.stickiness};
#else
  PriorityQueue pq;
#endif

  std::clog << "Using priority queue: " << pq.description() << '\n';

  thread_data = new ThreadData[settings.num_threads];
  std::seed_seq seq{settings.seed};
  auto seeds = new std::uint32_t[settings.num_threads];
  seq.generate(seeds, seeds + settings.num_threads);
  for (std::size_t i = 0; i < settings.num_threads; ++i) {
    thread_data[i].seed = seeds[i];
    thread_data[i].inserter =
        InsertingStrategy<key_type>(settings.insert_config, seeds[i] + 1);
  }
  delete[] seeds;
  thread_coordination::ThreadCoordinator coordinator{settings.num_threads};
  std::chrono::milliseconds duration;
  coordinator.run_task<Task>(std::ref(pq), std::ref(duration));
  coordinator.join();

  OperationCount sum_ops = std::transform_reduce(
      thread_data, thread_data + settings.num_threads, OperationCount{},
      std::plus<>(), [](auto const& d) { return d.op_count; });
  std::clog << "Insertions: " << sum_ops.num_insertions << '\n'
            << "Deletions: " << sum_ops.num_deletions << '\n'
            << "Failed deletions: " << sum_ops.num_failed_deletions << '\n';

#ifdef QUALITY
  std::cout << settings.num_threads << '\n';
  for (unsigned int t = 0; t < settings.num_threads; ++t) {
    for (auto const& [tick, key] : thread_data[t].ins_log) {
      std::cout << "i " << t << ' ' << tick << ' ' << key << '\n';
    }
  }

  for (unsigned int t = 0; t < settings.num_threads; ++t) {
    for (auto const& [tick, value] : thread_data[t].del_log) {
      std::cout << "d " << t << ' ' << tick << ' ' << get_thread_id(value)
                << ' ' << get_elem_id(value) << '\n';
    }
  }

  std::cout << std::flush;
#else
  std::cout << "Time (ms): " << std::chrono::milliseconds{duration}.count()
            << '\n'
            << "Ops/s: " << std::fixed << std::setprecision(1) << '\n'
            << static_cast<double>(settings.num_operations) /
                   std::chrono::duration<double>(duration).count()
            << std::endl;
#endif
  delete[] thread_data;
  delete[] prefill_keys;
  delete[] keys;
  return 0;
}
