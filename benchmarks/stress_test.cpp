#include "util/build_info.hpp"
#include "util/thread_coordination.hpp"
#include "wrapper/selector.hpp"

#include "cxxopts.hpp"

#ifdef LOG_OPERATIONS
#include "util/operation_log.hpp"
#endif

#ifdef WITH_PAPI
#include <papi.h>
#include <pthread.h>
#endif

#include <atomic>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <new>
#include <numeric>
#include <random>
#include <type_traits>
#include <utility>
#include <vector>

using key_type = unsigned long;
using value_type = unsigned long;

using pq_type = PQ<true, key_type, std::pair<key_type, value_type>>;
using handle_type = pq_type::handle_type;

struct Settings {
    enum class Mode { Update, Random, PopAll, PushPop };

    int num_threads = 4;
    long long prefill_per_thread = 1 << 20;
    long long iterations_per_thread = 1 << 22;
    Mode mode = Mode::Update;
    key_type min_prefill = 1;
    key_type max_prefill = 1 << 20;
    long min_update = 0;
    long max_update = 1 << 20;
    key_type min_random = 1;
    key_type max_random = 1 << 20;
    long long batch_size = 1 << 12;
    int seed = 1;
    std::chrono::seconds timeout{0};
#ifdef LOG_OPERATIONS
    std::filesystem::path log_file{"operation_log.txt"};
#endif
#ifdef WITH_PAPI
    std::vector<std::string> papi_events;
#endif

    static bool validate(Settings const& settings) {
        if (settings.num_threads <= 0) {
            std::cerr << "Error: Number of threads must be greater than 0\n";
            return false;
        }
        if (settings.min_prefill <= 0 || settings.min_random <= 0) {
            std::cerr << "Error: Keys must be greater than 0\n";
            return false;
        }
        if (settings.max_prefill < settings.min_prefill || settings.max_random < settings.min_random ||
            settings.max_update < settings.min_update) {
            std::cerr << "Error: Max must be greater than min\n";
            return false;
        }
        if (settings.batch_size <= 0) {
            std::cerr << "Error: batch size must be greater than 0\n";
            return false;
        }
        if (settings.prefill_per_thread < 0 || settings.iterations_per_thread < 0) {
            std::cerr << "Error: Prefill and iterations must be nonnegative\n";
            return false;
        }
#ifdef WITH_PAPI
        if (!settings.papi_events.empty()) {
            if (int ret = PAPI_library_init(PAPI_VER_CURRENT); ret != PAPI_VER_CURRENT) {
                std::cerr << "Error: Failed to initialize PAPI library\n";
                return false;
            }
            if (int ret = PAPI_thread_init(pthread_self); ret != PAPI_OK) {
                std::cerr << "Error: Failed to initialize PAPI thread support\n";
                return false;
            }
            for (auto const& name : settings.papi_events) {
                if (PAPI_query_named_event(name.c_str()) != PAPI_OK) {
                    std::cerr << "Error: PAPI event '" << name << "' not available\n";
                    return false;
                }
            }
        }
#endif
        return true;
    }
};

struct ThreadData {
    std::chrono::high_resolution_clock::time_point start_time{};
    std::chrono::high_resolution_clock::time_point end_time{};
    long long push_count = 0;
    long long pop_count = 0;
    long long failed_pop_count = 0;
#ifdef WITH_PAPI
    std::vector<long long> papi_counter{};
    bool papi_success = true;
#endif
#ifdef LOG_OPERATIONS
    struct PushLog {
        std::chrono::high_resolution_clock::time_point tick;
        std::pair<key_type, value_type> element;
    };
    struct PopLog {
        std::chrono::high_resolution_clock::time_point tick;
        value_type val;
    };
    std::vector<PushLog> pushes;
    std::vector<PopLog> pops;
#endif
};

class Context : public thread_coordination::Context {
    handle_type handle_;
    ThreadData data_;

   public:
    explicit Context(thread_coordination::Context ctx, handle_type handle)
        : thread_coordination::Context{std::move(ctx)}, handle_{std::move(handle)} {
    }

#ifdef LOG_OPERATIONS
    void push(std::pair<key_type, value_type> const& e) {
        handle_.push(e);
        auto tick = std::chrono::high_resolution_clock::now();
        data_.pushes.push_back({tick, e});
    }

    auto try_pop() {
        auto tick = std::chrono::high_resolution_clock::now();
        auto retval = handle_.try_pop();
        if (retval) {
            data_.pops.push_back({tick, retval->second});
        }
        return retval;
    }
#else
    void push(std::pair<key_type, value_type> e) {
        handle_.push(std::move(e));
    }

    auto try_pop() {
        return handle_.try_pop();
    }
#endif

    ThreadData& data() {
        return data_;
    }
};

class Update {
    std::vector<long long> data_;
    std::atomic_llong counter_{0};

   public:
    explicit Update(Settings const& settings)
        : data_(static_cast<std::size_t>(settings.iterations_per_thread * settings.num_threads)) {
    }

    template <typename RNG>
    void prepare(Settings const& settings, Context const& context, RNG&& rng) {
        std::generate_n(data_.begin() + context.id() * settings.iterations_per_thread, settings.iterations_per_thread,
                        [&rng, min = settings.min_update, max = settings.max_update]() {
                            return std::uniform_int_distribution<long>(min, max)(rng);
                        });
    }

    void work(Settings const& settings, Context& context) {
        auto offset = static_cast<value_type>(settings.num_threads * settings.prefill_per_thread);
        long long max = settings.iterations_per_thread * settings.num_threads;
        for (auto from = counter_.fetch_add(settings.batch_size, std::memory_order_relaxed); from < max;
             from = counter_.fetch_add(settings.batch_size, std::memory_order_relaxed)) {
            auto to = std::min(from + settings.batch_size, max);
            for (auto i = from; i < to; ++i) {
                auto key = context.try_pop();
                while (!key) {
                    ++context.data().failed_pop_count;
                    key = context.try_pop();
                }
                context.push(
                    {static_cast<key_type>(static_cast<long long>(key->first) + data_[static_cast<std::size_t>(i)]),
                     offset + static_cast<value_type>(i)});
            }
            context.data().push_count += to - from;
            context.data().pop_count += to - from;
        }
    }

    std::size_t required_capacity(Settings const& settings) const noexcept {
        return static_cast<std::size_t>(settings.num_threads);
    }

    std::size_t pushes_per_thread(Settings const& settings) const noexcept {
        return static_cast<std::size_t>(settings.iterations_per_thread);
    }

    std::size_t pops_per_thread(Settings const& settings) const noexcept {
        return static_cast<std::size_t>(settings.iterations_per_thread);
    }
};

class Random {
    std::vector<key_type> keys_;
    std::atomic_llong counter_{0};

   public:
    explicit Random(Settings const& settings)
        : keys_(static_cast<std::size_t>(settings.iterations_per_thread * settings.num_threads)) {
    }

    template <typename RNG>
    void prepare(Settings const& settings, Context const& context, RNG&& rng) {
        std::generate_n(keys_.begin() + context.id() * settings.iterations_per_thread, settings.iterations_per_thread,
                        [&rng, min = settings.min_random, max = settings.max_random]() {
                            return std::uniform_int_distribution<key_type>(min, max)(rng);
                        });
    }

    void work(Settings const& settings, Context& context) {
        auto offset = static_cast<value_type>(settings.num_threads * settings.prefill_per_thread);
        long long max = settings.iterations_per_thread * settings.num_threads;
        for (auto from = counter_.fetch_add(settings.batch_size, std::memory_order_relaxed); from < max;
             from = counter_.fetch_add(settings.batch_size, std::memory_order_relaxed)) {
            auto to = std::min(from + settings.batch_size, max);
            for (auto i = from; i < to; ++i) {
                while (!context.try_pop()) {
                    ++context.data().failed_pop_count;
                }
                context.push({keys_[static_cast<std::size_t>(i)], offset + static_cast<value_type>(i)});
            }
            context.data().push_count += to - from;
            context.data().pop_count += to - from;
        }
    }

    std::size_t required_capacity(Settings const& settings) const noexcept {
        return static_cast<std::size_t>(settings.num_threads);
    }

    std::size_t pushes_per_thread(Settings const& settings) const noexcept {
        return static_cast<std::size_t>(settings.iterations_per_thread);
    }

    std::size_t pops_per_thread(Settings const& settings) const noexcept {
        return static_cast<std::size_t>(settings.iterations_per_thread);
    }
};

class PopAll {
    std::atomic_llong counter_{0};

   public:
    template <typename RNG>
    void prepare(Settings const&, Context const&, RNG const&) {
    }

    void work(Settings const& settings, Context& context) {
        while (true) {
            long long count = 0;
            while (context.try_pop()) {
                ++count;
            }
            context.data().pop_count += count;
            auto current = counter_.fetch_add(count, std::memory_order_relaxed) + count;
            if (current >= settings.iterations_per_thread * settings.num_threads) {
                break;
            }
        }
    }

    std::size_t required_capacity(Settings const&) const noexcept {
        return 0;
    }

    std::size_t pushes_per_thread(Settings const&) const noexcept {
        return 0;
    }

    std::size_t pops_per_thread(Settings const& settings) const noexcept {
        return static_cast<std::size_t>(settings.prefill_per_thread);
    }
};

class PushPop {
    std::vector<key_type> keys_;
    std::atomic_llong push_counter_{0};
    std::atomic_llong pop_counter_{0};

   public:
    explicit PushPop(Settings const& settings)
        : keys_(static_cast<std::size_t>(settings.iterations_per_thread * settings.num_threads)) {
    }

    template <typename RNG>
    void prepare(Settings const& settings, Context const& context, RNG&& rng) {
        std::generate_n(keys_.begin() + context.id() * settings.iterations_per_thread, settings.iterations_per_thread,
                        [&rng, min = settings.min_random, max = settings.max_random]() {
                            return std::uniform_int_distribution<key_type>(min, max)(rng);
                        });
    }

    void work(Settings const& settings, Context& context) {
        auto offset = static_cast<value_type>(settings.num_threads * settings.prefill_per_thread);
        long long max = settings.iterations_per_thread * settings.num_threads;
        for (auto from = push_counter_.fetch_add(settings.batch_size, std::memory_order_relaxed); from < max;
             from = push_counter_.fetch_add(settings.batch_size, std::memory_order_relaxed)) {
            auto to = std::min(from + settings.batch_size, max);
            for (auto i = from; i < to; ++i) {
                context.push({keys_[static_cast<std::size_t>(i)], offset + static_cast<value_type>(i)});
            }
            context.data().push_count += to - from;
        }
        for (auto from = pop_counter_.fetch_add(settings.batch_size, std::memory_order_relaxed); from < max;
             from = pop_counter_.fetch_add(settings.batch_size, std::memory_order_relaxed)) {
            auto to = std::min(from + settings.batch_size, max);
            for (auto i = from; i < to; ++i) {
                while (!context.try_pop()) {
                    ++context.data().failed_pop_count;
                }
            }
            context.data().pop_count += to - from;
        }
    }

    std::size_t required_capacity(Settings const& settings) const noexcept {
        return static_cast<std::size_t>(settings.iterations_per_thread * settings.num_threads);
    }

    std::size_t pushes_per_thread(Settings const& settings) const noexcept {
        return static_cast<std::size_t>(settings.iterations_per_thread);
    }

    std::size_t pops_per_threads(Settings const& settings) const noexcept {
        return static_cast<std::size_t>(settings.iterations_per_thread);
    }
};

static constexpr auto mode_name(Settings::Mode mode) noexcept {
    switch (mode) {
        case Settings::Mode::Update:
            return "update";
        case Settings::Mode::Random:
            return "random";
        case Settings::Mode::PopAll:
            return "pop_all";
        case Settings::Mode::PushPop:
            return "push_pop";
        default:
            return "unknown";
    }
};

void output_json(pq_type const& pq, Settings const& settings, std::vector<ThreadData> const& data) {
    auto print_kv = [](auto const& key, auto const& v) { std::cout << std::quoted(key) << ':' << v; };
    auto last_end = std::max_element(data.begin(), data.end(), [](auto const& lhs, auto const& rhs) {
                        return lhs.end_time < rhs.end_time;
                    })->end_time;
    auto first_start = std::min_element(data.begin(), data.end(), [](auto const& lhs, auto const& rhs) {
                           return lhs.start_time < rhs.start_time;
                       })->start_time;
    std::cout << '{';
    print_kv("pq", pq.describe_json());
    std::cout << ',' << std::quoted("settings") << ':' << '{';
    print_kv("num_threads", settings.num_threads);
    std::cout << ',';
    print_kv("prefill_per_thread", settings.prefill_per_thread);
    std::cout << ',';
    print_kv("iterations_per_thread", settings.iterations_per_thread);
    std::cout << ',';
    print_kv("prefill_min", settings.min_prefill);
    std::cout << ',';
    print_kv("prefill_max", settings.max_prefill);
    std::cout << ',';
    print_kv("mode", std::quoted(mode_name(settings.mode)));
    if (settings.mode == Settings::Mode::Update) {
        std::cout << ',';
        print_kv("update_min", settings.min_update);
        std::cout << ',';
        print_kv("update_max", settings.max_update);
    } else if (settings.mode == Settings::Mode::Random) {
        std::cout << ',';
        print_kv("random_min", settings.min_random);
        std::cout << ',';
        print_kv("random_max", settings.max_random);
    }
    std::cout << ',';
    print_kv("batch_size", settings.batch_size);
    std::cout << ',';
    print_kv("seed", settings.seed);
    std::cout << ',';
    print_kv("timeout", settings.timeout.count());
#ifdef WITH_PAPI
    std::cout << ',';
    std::cout << std::quoted("papi_events") << ':' << '[';
    for (std::size_t i = 0; i < settings.papi_events.size(); ++i) {
        std::cout << std::quoted(settings.papi_events[i]);
        if (i != settings.papi_events.size() - 1) {
            std::cout << ',';
        }
    }
    std::cout << ']';
#endif
    std::cout << '}';
    std::cout << ',' << std::quoted("results") << ':' << '{';
    print_kv("total_time", (last_end - first_start).count());
    std::cout << ',' << std::quoted("thread_stats") << ':' << '[';
    for (auto it = data.begin(); it != data.end(); ++it) {
        std::cout << '{';
        print_kv("time", (it->end_time - it->start_time).count());
        std::cout << ',';
        print_kv("pushes", it->push_count);
        std::cout << ',';
        print_kv("pops", it->pop_count);
        std::cout << ',';
        print_kv("failed_pops", it->failed_pop_count);
#ifdef WITH_PAPI
        std::cout << ',' << std::quoted("papi_events") << ':';
        if (it->papi_success) {
            std::cout << '[';
            for (std::size_t i = 0; i < settings.papi_events.size(); ++i) {
                std::cout << '{';
                print_kv("name", std::quoted(settings.papi_events[i]));
                std::cout << ',';
                print_kv("count", it->papi_counter[i]);
                std::cout << '}';
                if (i != settings.papi_events.size() - 1) {
                    std::cout << ',';
                }
            }
            std::cout << ']';
        } else {
            std::cout << "null";
        }
#endif
        std::cout << '}';
        if (it != std::prev(data.end())) {
            std::cout << ',';
        }
    }
    std::cout << ']' << '}' << '}' << std::endl;
}

#ifdef WITH_PAPI
int prepare_papi(Settings const& settings) {
    if (int ret = PAPI_register_thread(); ret != PAPI_OK) {
        throw std::runtime_error{"Failed to register thread for PAPI"};
    }
    int event_set = PAPI_NULL;
    if (int ret = PAPI_create_eventset(&event_set); ret != PAPI_OK) {
        throw std::runtime_error{"Failed to create PAPI event set"};
    }
    for (auto const& name : settings.papi_events) {
        auto event = PAPI_NULL;
        if (PAPI_event_name_to_code(name.c_str(), &event) != PAPI_OK) {
            throw std::runtime_error{"Failed to resolve PAPI event '" + name + '\''};
        }
        if (PAPI_add_event(event_set, event) != PAPI_OK) {
            throw std::runtime_error{"Failed to add PAPI event '" + name + '\''};
        }
    }
    return event_set;
}
#endif

template <typename Task>
ThreadData benchmark_thread(thread_coordination::Context thread_context, Settings const& settings, pq_type& pq,
                            Task& task) {
    Context context{std::move(thread_context), pq.get_handle()};
#ifdef WITH_PAPI
    int event_set = PAPI_NULL;
    if (!settings.papi_events.empty()) {
        context.data().papi_counter.resize(settings.papi_events.size());
        try {
            event_set = prepare_papi(settings);
        } catch (std::exception const& e) {
            context.write(std::cerr) << e.what() << '\n';
        }
    }
#endif
#ifdef LOG_OPERATIONS
    context.data().pushes.reserve(static_cast<std::size_t>(settings.prefill_per_thread) +
                                  2 * task.pushes_per_thread(settings));
    data.pops.reserve(2 * task.pops_per_thread(settings));
#endif

    {
        std::vector<key_type> prefill(static_cast<std::size_t>(settings.prefill_per_thread));

        if (context.id() == 0) {
            std::clog << "Preparing...\n";
        }
        {
            std::seed_seq seed{settings.seed, context.id()};
            std::default_random_engine rng(seed);
            context.synchronize();
            std::generate(prefill.begin(), prefill.end(),
                          [&rng, min = settings.min_prefill, max = settings.max_prefill]() {
                              return std::uniform_int_distribution<key_type>(min, max)(rng);
                          });
            context.synchronize();
            task.prepare(settings, context, rng);
        }
        context.synchronize();
        if (context.id() == 0) {
            std::clog << "Prefilling...\n";
        }
        context.synchronize();
        for (auto i = 0LL; i < settings.prefill_per_thread; ++i) {
            context.push({prefill[static_cast<std::size_t>(i)],
                          static_cast<value_type>(context.id() * settings.prefill_per_thread + i)});
        }
    }
    context.synchronize();
    if (context.id() == 0) {
        std::clog << "Running benchmark...\n";
    }
    context.synchronize();
#ifdef WITH_PAPI
    if (!settings.papi_events.empty()) {
        if (int ret = PAPI_start(event_set); ret != PAPI_OK) {
            context.write(std::cerr) << "Failed to start performance counters\n";
        }
    }
#endif
    context.data().start_time = std::chrono::high_resolution_clock::now();
    task.work(settings, context);
    context.data().end_time = std::chrono::high_resolution_clock::now();
#ifdef WITH_PAPI
    if (!settings.papi_events.empty()) {
        if (int ret = PAPI_stop(event_set, context.data().papi_counter.data()); ret != PAPI_OK) {
            context.write(std::cerr) << "Failed to stop performance counters\n";
        }
    }
#endif
    return std::move(context.data());
}

#ifdef LOG_OPERATIONS
void output_log(std::ostream& log_out, std::vector<ThreadData> const& thread_data) {
    std::vector<ThreadData::PushLog> pushes;
    pushes.reserve(std::accumulate(thread_data.begin(), thread_data.end(), 0UL,
                                   [](std::size_t sum, auto const& e) { return sum + e.pushes.size(); }));
    std::vector<ThreadData::PopLog> pops;
    pushes.reserve(std::accumulate(thread_data.begin(), thread_data.end(), 0UL,
                                   [](std::size_t sum, auto const& e) { return sum + e.pops.size(); }));
    for (auto const& e : thread_data) {
        pushes.insert(pushes.end(), e.pushes.begin(), e.pushes.end());
        pops.insert(pops.end(), e.pops.begin(), e.pops.end());
    }
    std::sort(pushes.begin(), pushes.end(), [](auto const& lhs, auto const& rhs) { return lhs.tick < rhs.tick; });
    std::vector<std::size_t> push_index(pushes.size());
    for (std::size_t i = 0; i < pushes.size(); ++i) {
        push_index[static_cast<std::size_t>(pushes[i].element.second)] = i;
    }
    std::sort(pops.begin(), pops.end(), [](auto const& lhs, auto const& rhs) { return lhs.tick < rhs.tick; });
    std::size_t i = 0;
    for (auto const& pop : pops) {
        auto ref_index = push_index[static_cast<std::size_t>(pop.val)];
        while ((i != pushes.size() && pushes[i].tick < pop.tick) || i <= ref_index) {
            log_out << i << ' ' << pushes[i].element.first << '\n';
            ++i;
        }
        log_out << '-' << ref_index << '\n';
    }
    for (; i < pushes.size(); ++i) {
        log_out << i << ' ' << pushes[i].element.first << '\n';
    }
}
#endif

bool run_benchmark(Settings const& settings, cxxopts::ParseResult const& parse_result) {
    if (!Settings::validate(settings)) {
        return false;
    }

#ifdef LOG_OPERATIONS
    std::ofstream log_out;
    if (!settings.log_file.empty()) {
        log_out = std::ofstream(settings.log_file);
    }
    if (!log_out) {
        std::cerr << "Error: Could not open file " << settings.log_file << " for writing" << std::endl;
        return false;
    }
#endif

    std::vector<ThreadData> thread_data(static_cast<std::size_t>(settings.num_threads));
    std::unique_ptr<pq_type> pq;
    auto dispatch = [&](auto& task) {
        pq = std::make_unique<pq_type>(settings.num_threads,
                                       static_cast<std::size_t>(settings.prefill_per_thread * settings.num_threads) +
                                           task.required_capacity(settings),
                                       parse_result);
        auto dispatcher = thread_coordination::Dispatcher(
            affinity::NUMA{cores_per_numa_node, num_numa_nodes}, settings.num_threads, [&](auto t_ctx) {
                thread_data[static_cast<std::size_t>(t_ctx.id())] =
                    benchmark_thread(std::move(t_ctx), settings, *pq, task);
            });
        dispatcher.wait();
    };

    switch (settings.mode) {
        case Settings::Mode::Update: {
            auto task = Update(settings);
            dispatch(task);
            break;
        }
        case Settings::Mode::Random: {
            auto task = Random(settings);
            dispatch(task);
            break;
        }
        case Settings::Mode::PopAll: {
            auto task = PopAll();
            dispatch(task);
            break;
        }
        case Settings::Mode::PushPop: {
            auto task = PushPop(settings);
            dispatch(task);
            break;
        }
    }

    output_json(*pq, settings, thread_data);
#ifdef LOG_OPERATIONS
    output_log(log_out, thread_data);
#endif
    std::clog << "Finished\n";
    return true;
}

int main(int argc, char* argv[]) {
    write_build_info(std::clog);
    std::clog << "\nCommand line:";
    for (int i = 0; i < argc; ++i) {
        std::clog << ' ' << argv[i];
    }
    std::clog << '\n';

    Settings settings;
    cxxopts::Options cmd(argv[0]);
    int timeout_ms = 0;
    cmd.add_options()
        // clang-format off
        ("h,help", "Print this help", cxxopts::value<bool>())
        ("j,threads", "Number of threads", cxxopts::value<int>(settings.num_threads), "NUMBER")
        ("p,prefill", "Prefill per thread", cxxopts::value<long long>(settings.prefill_per_thread), "NUMBER")
        ("n,iterations", "Number of iterations", cxxopts::value<long long>(settings.iterations_per_thread), "NUMBER")
        ("mode", "Operation mode ([u]pdate, [r]andom, [p]op, pu[s]h, push [a]scending", cxxopts::value<char>(), "STRING")
        ("min-prefill", "Min prefill key", cxxopts::value<key_type>(settings.min_prefill), "NUMBER")
        ("max-prefill", "Max prefill key", cxxopts::value<key_type>(settings.max_prefill), "NUMBER")
        ("min-update", "Min update", cxxopts::value<long>(settings.min_update), "NUMBER")
        ("max-update", "Max update", cxxopts::value<long>(settings.max_update), "NUMBER")
        ("min-random", "Min random key", cxxopts::value<key_type>(settings.min_random), "NUMBER")
        ("max-random", "Max random key", cxxopts::value<key_type>(settings.max_random), "NUMBER")
        ("batch-size", "Batch size", cxxopts::value<long long>(settings.batch_size), "NUMBER")
        ("s,seed", "Initial seed", cxxopts::value<int>(settings.seed), "NUMBER")
        ("timeout", "Timeout in milliseconds", cxxopts::value<int>(timeout_ms), "NUMBER")
#ifdef LOG_OPERATIONS
        ("log-file", "File to write the operation log to", cxxopts::value<std::filesystem::path>(settings.log_file), "PATH")
#endif
#ifdef WITH_PAPI
        ("r,pc", "Performance counters", cxxopts::value<std::vector<std::string>>(settings.papi_events))
#endif
        // clang-format on
        ;
    pq_type::add_options(cmd);
#ifdef LOG_OPERATIONS
    cmd.parse_positional({"log-file"});
#endif
    cxxopts::ParseResult args;
    try {
        args = cmd.parse(argc, argv);
    } catch (std::exception const& e) {
        std::cerr << "Error parsing command line: " << e.what() << std::endl;
        std::cerr << cmd.help() << std::endl;
        return EXIT_FAILURE;
    }

    if (args.count("help") > 0) {
        std::clog << cmd.help() << std::endl;
        return EXIT_SUCCESS;
    }
    if (args.count("mode") > 0) {
        switch (args["mode"].as<char>()) {
            case 'u':
                settings.mode = Settings::Mode::Update;
                break;
            case 'r':
                settings.mode = Settings::Mode::Random;
                break;
            case 'p':
                settings.mode = Settings::Mode::PopAll;
                break;
            case 'a':
                settings.mode = Settings::Mode::PushPop;
                break;
            default:
                std::cerr << "Error parsing command line: Invalid mode\n";
                std::cerr << cmd.help() << std::endl;
                return EXIT_FAILURE;
        }
    }
    settings.timeout = std::chrono::seconds(timeout_ms);

    bool success = run_benchmark(settings, args);
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
