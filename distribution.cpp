#include "cxxopts.hpp"
#include "multiqueue/configurations.hpp"
#include "multiqueue/int_multiqueue.hpp"
#include "system_config.hpp"
#include "utils/inserting_strategy.hpp"
#include "utils/thread_coordination.hpp"
#include "utils/threading.hpp"

#include <x86intrin.h>
#include <array>
#include <atomic>
#include <cassert>
#include <chrono>
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

using key_type = unsigned long;
using value_type = unsigned long;
static_assert(std::is_unsigned_v<value_type>, "Value type must be unsigned");

using namespace std::chrono_literals;

struct Settings {
  std::size_t prefill_size = 1'000'000;
  std::chrono::nanoseconds sleep_between_operations = 0ns;
  std::size_t num_operations_between_samples = 100'000;
  std::size_t num_samples = 10;
  std::size_t sample_size = 0;
  unsigned int num_threads = 4;
  InsertConfig<key_type> insert_config{
      InsertPolicy::Uniform,
      KeyDistribution::Uniform,
      std::numeric_limits<value_type>::min(),
      std::numeric_limits<value_type>::max() - 3,  // Some pqs use sentinels
      1,
      100,
  };
  std::uint32_t seed;
};

std::uint32_t* thread_seeds;

std::atomic_bool start_flag;

std::vector<std::vector<std::size_t>> samples;
std::atomic_uint64_t operation_counter;

template <typename PriorityQueue>
struct Task {
  static void run(thread_coordination::Context ctx, PriorityQueue& pq,
                  Settings const& settings) {
    unsigned int stage = 0;

    std::seed_seq seq{thread_seeds[ctx.get_id()]};
    auto gen = std::mt19937(seq);
    auto dist = std::uniform_int_distribution<long>(
        100, settings.sleep_between_operations.count());
    auto handle = pq.get_handle(ctx.get_id());

    auto inserter = InsertingStrategy<key_type>{
        ctx.get_id(), settings.insert_config, thread_seeds[ctx.get_id()] + 1};

    if (ctx.is_main()) {
      if (settings.prefill_size > 0) {
        std::clog << "Prefilling..." << std::flush;
        for (size_t i = 0; i < settings.prefill_size; ++i) {
          key_type const key = inserter.get_key();
          pq.push(handle, {key, key});
        }
        std::clog << "done" << std::endl;
      }
      std::clog << "Starting test...\n";
      std::clog << "Taking sample 1" << std::endl;
      if (settings.sample_size == 0) {
        samples.push_back(pq.get_distribution());
      } else {
        samples.push_back(pq.get_top_distribution(settings.sample_size));
      }
    }
    ctx.synchronize(stage++, [&ctx]() { ctx.notify_coordinator(); });
    while (!start_flag.load(std::memory_order_relaxed)) {
      _mm_pause();
    }
    std::atomic_thread_fence(std::memory_order_acquire);
    std::pair<key_type, value_type> retval;
    for (std::size_t i = 1; i < settings.num_samples; ++i) {
      ctx.synchronize(stage++);
      while (operation_counter.load(std::memory_order_relaxed) <
             settings.num_operations_between_samples) {
        if (inserter.insert()) {
          key_type const key = inserter.get_key();
          pq.push(handle, {key, key});
          operation_counter.fetch_add(1, std::memory_order_relaxed);
        } else {
          if (pq.extract_top(handle, retval)) {
            operation_counter.fetch_add(1, std::memory_order_relaxed);
          }
        }
        if (settings.sleep_between_operations > 0us) {
          std::this_thread::sleep_for(std::chrono::nanoseconds{dist(gen)});
        }
      }
      ctx.synchronize(stage++, [&pq, &settings]() {
        std::clog << "Taking sample " << samples.size() + 1 << std::endl;
        if (settings.sample_size == 0) {
          samples.push_back(pq.get_distribution());
        } else {
          samples.push_back(pq.get_top_distribution(settings.sample_size));
        }
        operation_counter = 0;
      });
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

template <typename PriorityQueue>
void benchmark(Settings const& settings) {
  std::clog << "Using priority queue: " << PriorityQueue::description() << '\n';
  PriorityQueue pq{settings.num_threads, settings.seed};
  thread_coordination::ThreadCoordinator coordinator{settings.num_threads};
  coordinator.run<Task<PriorityQueue>>(std::ref(pq), settings);
  coordinator.wait_until_notified();
  start_flag.store(true, std::memory_order_release);
  coordinator.join();
}

template <unsigned int c, unsigned int k>
struct MQConfiguration : multiqueue::configuration::NoBuffering {
  static constexpr unsigned int C = c;
  static constexpr unsigned int K = k;
};

template <unsigned int c>
bool dispatch_k(unsigned int k, Settings const& settings) {
  switch (k) {
    case 1:
      benchmark<multiqueue::int_multiqueue<key_type, value_type,
                                           MQConfiguration<c, 1>>>(settings);
      return true;
    case 4:
      benchmark<multiqueue::int_multiqueue<key_type, value_type,
                                           MQConfiguration<c, 4>>>(settings);
      return true;
    case 8:
      benchmark<multiqueue::int_multiqueue<key_type, value_type,
                                           MQConfiguration<c, 8>>>(settings);
      return true;
    case 16:
      benchmark<multiqueue::int_multiqueue<key_type, value_type,
                                           MQConfiguration<c, 16>>>(settings);
      return true;
    default:
      std::cerr << "Invalid stickiness!\n";
  }
  return false;
}

bool dispatch_c(unsigned int c, unsigned int k, Settings const& settings) {
  switch (c) {
    case 2:
      return dispatch_k<2>(k, settings);
    case 4:
      return dispatch_k<4>(k, settings);
    case 8:
      return dispatch_k<8>(k, settings);
    case 16:
      return dispatch_k<16>(k, settings);
    default:
      std::cerr << "Invalid number of queues per thread!\n";
  }
  return false;
}

int main(int argc, char* argv[]) {
  Settings settings{};

  cxxopts::Options options(
      "Distribution",
      "This executable records the distribution of keys in the multiqueue");
  // clang-format off
    options.add_options()
      ("c,queues", "The number of queues per thread"
       "(default: 4)", cxxopts::value<unsigned int>(), "NUMBER")
      ("k,stickiness", "The stickiness parameter"
       "(default: 1)", cxxopts::value<unsigned int>(), "NUMBER")
      ("n,prefill", "Specify the number of elements to prefill the queue with "
       "(default: 1'000'000)", cxxopts::value<size_t>(), "NUMBER")
      ("i,insert", "Specify the insert policy as one of \"uniform\", \"split\", \"producer\", \"alternating\" "
       "(default: uniform)", cxxopts::value<std::string>(), "ARG")
      ("j,threads", "Specify the number of threads "
       "(default: 4)", cxxopts::value<unsigned int>(), "NUMBER")
      ("w,sleep", "Specify the sleep time between operations in ns"
       "(default: 0)", cxxopts::value<unsigned int>(), "NUMBER")
      ("d,distribution", "Specify the key distribution as one of \"uniform\", \"dijkstra\", \"ascending\", \"descending\", \"threadid\" "
       "(default: uniform)", cxxopts::value<std::string>(), "ARG")
      ("m,max", "Specify the max key "
       "(default: MAX)", cxxopts::value<key_type>(), "NUMBER")
      ("l,min", "Specify the min key "
       "(default: 0)", cxxopts::value<key_type>(), "NUMBER")
      ("s,seed", "Specify the initial seed"
       "(default: 0)", cxxopts::value<std::uint32_t>(), "NUMBER")
      ("o,operations", "Specify the number of operations between samples"
       "(default: 10'000)", cxxopts::value<std::size_t>(), "NUMBER")
      ("p,samples", "Specify the number of samples"
       "(default: 10)", cxxopts::value<std::size_t>(), "NUMBER")
      ("z,samplesize", "Specify the number of keys to sample. 0 for all elements in the queue"
       "(default: 0)", cxxopts::value<std::size_t>(), "NUMBER")
      ("h,help", "Print this help");
  // clang-format on

  unsigned int c = 4;
  unsigned int k = 1;
  try {
    auto result = options.parse(argc, argv);
    if (result.count("help")) {
      std::cerr << options.help() << std::endl;
      return 0;
    }
    if (result.count("prefill") > 0) {
      settings.prefill_size = result["prefill"].as<size_t>();
    }
    if (result.count("insert") > 0) {
      std::string policy = result["insert"].as<std::string>();
      if (policy == "uniform") {
        settings.insert_config.insert_policy = InsertPolicy::Uniform;
      } else if (policy == "split") {
        settings.insert_config.insert_policy = InsertPolicy::Split;
      } else if (policy == "producer") {
        settings.insert_config.insert_policy = InsertPolicy::Producer;
      } else if (policy == "alternating") {
        settings.insert_config.insert_policy = InsertPolicy::Alternating;
      } else {
        std::cerr << "Unknown insert policy \"" << policy << "\"\n";
        return 1;
      }
    }
    if (result.count("threads") > 0) {
      settings.num_threads = result["threads"].as<unsigned int>();
    }
    if (result.count("sleep") > 0) {
      settings.sleep_between_operations =
          std::chrono::nanoseconds{result["sleep"].as<unsigned int>()};
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
      } else if (dist == "threadid") {
        settings.insert_config.key_distribution = KeyDistribution::ThreadId;
      } else {
        std::cerr << "Unknown key distribution \"" << dist << "\"\n";
        return 1;
      }
    }
    if (result.count("max") > 0) {
      settings.insert_config.max_key = result["max"].as<key_type>();
    }
    if (result.count("min") > 0) {
      settings.insert_config.min_key = result["min"].as<key_type>();
    }
    if (result.count("operations") > 0) {
      settings.num_operations_between_samples =
          result["operations"].as<std::size_t>();
    }
    if (result.count("samples") > 0) {
      settings.num_samples = result["samples"].as<std::size_t>();
    }
    if (result.count("samplesize") > 0) {
      settings.sample_size = result["samplesize"].as<std::size_t>();
    }
    if (result.count("seed") > 0) {
      settings.seed = result["seed"].as<std::uint32_t>();
    }
    if (result.count("queues") > 0) {
      c = result["queues"].as<unsigned int>();
    }
    if (result.count("stickiness") > 0) {
      k = result["stickiness"].as<unsigned int>();
    }
  } catch (cxxopts::OptionParseException const& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }

#ifndef NDEBUG
  std::clog << "Using debug build!\n\n";
#endif

  std::clog << "Settings: \n\t"
            << "Prefill size: " << settings.prefill_size << "\n\t"
            << "Sleep between operations: "
            << settings.sleep_between_operations.count() << " ns\n\t"
            << "Num operations between samples: "
            << settings.num_operations_between_samples << "\n\t"
            << "Num samples: " << settings.num_samples << "\n\t"
            << "Samplesize: " << settings.sample_size << "\n\t"
            << "Threads: " << settings.num_threads << "\n\t"
            << "Insert policy: "
            << get_insert_policy_name(settings.insert_config.insert_policy)
            << "\n\t"
            << "Min key: " << settings.insert_config.min_key << "\n\t"
            << "Max key: " << settings.insert_config.max_key << "\n\t"
            << "Key distribution: "
            << get_key_distribution_name(
                   settings.insert_config.key_distribution)
            << "\n\t"
            << "Dijkstra min increase: "
            << settings.insert_config.dijkstra_min_increase << "\n\t"
            << "Dijkstra max increase: "
            << settings.insert_config.dijkstra_max_increase << "\n\t"
            << "Seed: " << settings.seed;
  std::clog << "\n\n";

  operation_counter = 0;

  std::seed_seq seq{settings.seed + 1};
  thread_seeds = new std::uint32_t[settings.num_threads];
  seq.generate(thread_seeds, thread_seeds + settings.num_threads);
  start_flag.store(false, std::memory_order_relaxed);
  std::atomic_thread_fence(std::memory_order_release);
  bool success = dispatch_c(c, k, settings);
  if (!success) {
    return 1;
  }

  std::cout << samples.size() << ' ' << samples[0].size() << '\n';
  for (auto const& s : samples) {
    std::copy(s.begin(), s.end() - 1,
              std::ostream_iterator<std::size_t>(std::cout, " "));
    std::cout << s.back();
    std::cout << '\n';
  }

  std::cout << std::flush;

  delete[] thread_seeds;
  return 0;
}
