#include "util/build_info.hpp"
#include "util/thread_coordination.hpp"
#include "wrapper/selector.hpp"
#ifdef LOG_OPERATIONS
#include "util/operation_log.hpp"
#endif

#include <cxxopts.hpp>

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
#ifdef LOG_OPERATIONS
    std::ofstream log_file;
#endif
#ifdef WITH_PAPI
    std::vector<std::string> papi_events;
#endif

    static void add_cmd_options(cxxopts::Options& cmd) {
        // clang-format off
    cmd.add_options()
        ("j,threads", "Number of threads", cxxopts::value<int>(), "NUMBER")
        ("p,prefill", "Prefill per thread", cxxopts::value<long long>(), "NUMBER")
        ("n,iterations", "Number of iterations per thread", cxxopts::value<long long>(), "NUMBER")
        ("mode", "Operation mode ([u]pdate, [r]andom, pop [a]ll, [p]ush pop)", cxxopts::value<char>(), "STRING")
        ("min-prefill", "Min prefill key", cxxopts::value<key_type>(), "NUMBER")
        ("max-prefill", "Max prefill key", cxxopts::value<key_type>(), "NUMBER")
        ("min-update", "Min update", cxxopts::value<long>(), "NUMBER")
        ("max-update", "Max update", cxxopts::value<long>(), "NUMBER")
        ("min-random", "Min random key", cxxopts::value<key_type>(), "NUMBER")
        ("max-random", "Max random key", cxxopts::value<key_type>(), "NUMBER")
        ("batch-size", "Batch size", cxxopts::value<long long>(), "NUMBER")
        ("s,seed", "Initial seed", cxxopts::value<int>(), "NUMBER")
#ifdef LOG_OPERATIONS
        ("log-file", "File to write the operation log to", cxxopts::value<std::filesystem::path>(), "PATH")
#endif
#ifdef WITH_PAPI
        ("r,pc", "Performance counters", cxxopts::value<std::vector<std::string>>())
#endif
        ;
        // clang-format on
#ifdef LOG_OPERATIONS
        cmd.parse_positional({"log-file"});
#endif
    }

    static std::optional<Settings> from_cmd(cxxopts::ParseResult const& args) {
        Settings settings;
        if (args.count("threads") > 0) {
            settings.num_threads = args["threads"].as<int>();
            if (settings.num_threads <= 0) {
                std::cerr << "Error: Number of threads must be greater than 0\n";
                return {};
            }
        }
        if (args.count("prefill") > 0) {
            settings.prefill_per_thread = args["prefill"].as<long long>();
            if (settings.prefill_per_thread < 0) {
                std::cerr << "Error: Prefill must be nonnegative\n";
                return {};
            }
        }
        if (args.count("iterations") > 0) {
            settings.iterations_per_thread = args["iterations"].as<long long>();
            if (settings.iterations_per_thread < 0) {
                std::cerr << "Error: Iterations must be nonnegative\n";
                return {};
            }
        }
        if (args.count("mode") > 0) {
            switch (args["mode"].as<char>()) {
                case 'u':
                    settings.mode = Mode::Update;
                    break;
                case 'r':
                    settings.mode = Mode::Random;
                    break;
                case 'a':
                    settings.mode = Mode::PopAll;
                    break;
                case 'p':
                    settings.mode = Mode::PushPop;
                    break;
                default:
                    std::cerr << "Error: Invalid mode\n";
                    return {};
            }
        }
        if (args.count("min-prefill") > 0) {
            settings.min_prefill = args["min-prefill"].as<key_type>();
            if (settings.min_prefill <= 0) {
                std::cerr << "Error: Prefill keys must be greater than 0\n";
                return {};
            }
        }
        if (args.count("max-prefill") > 0) {
            settings.max_prefill = args["max-prefill"].as<key_type>();
            if (settings.max_prefill < settings.min_prefill) {
                std::cerr << "Error: Invalid prefill range\n";
                return {};
            }
        }
        if (args.count("min-update") > 0) {
            settings.min_update = args["min-update"].as<long>();
        }
        if (args.count("max-update") > 0) {
            settings.max_update = args["max-update"].as<long>();
            if (settings.max_update < settings.min_update) {
                std::cerr << "Error: Invalid update range\n";
                return {};
            }
        }
        if (args.count("min-random") > 0) {
            settings.min_random = args["min-random"].as<key_type>();
            if (settings.min_random <= 0) {
                std::cerr << "Error: Random keys must be greater than 0\n";
                return {};
            }
        }
        if (args.count("max-random") > 0) {
            settings.max_random = args["max-random"].as<key_type>();
            if (settings.max_random < settings.min_random) {
                std::cerr << "Error: Invalid random range\n";
                return {};
            }
        }
        if (args.count("batch-size") > 0) {
            settings.batch_size = args["batch-size"].as<long long>();
            if (settings.batch_size <= 0) {
                std::cerr << "Error: batch size must be greater than 0\n";
                return {};
            }
        }
        if (args.count("seed") > 0) {
            settings.seed = args["seed"].as<int>();
        }
#ifdef LOG_OPERATIONS
        std::string log_file_name =
            args.count("log-file") > 0 ? args["log-file"].as<std::filesystem::path>() : "operation_log.txt";
        settings.log_out.open(log_file_name);
        if (!settings.log_out) {
            std::cerr << "Error: Could not open file '" << log_file_name << "' for writing" << std::endl;
            return {};
        }
#endif
#ifdef WITH_PAPI
        if (args.count("pc") > 0) {
            settings.papi_events = args["pc"].as<std::vector<std::string>>();
            if (!settings.papi_events.empty()) {
                if (int ret = PAPI_library_init(PAPI_VER_CURRENT); ret != PAPI_VER_CURRENT) {
                    std::cerr << "Error: Failed to initialize PAPI library\n";
                    return {};
                }
                if (int ret = PAPI_thread_init(pthread_self); ret != PAPI_OK) {
                    std::cerr << "Error: Failed to initialize PAPI thread support\n";
                    return {};
                }
                for (auto const& name : settings.papi_events) {
                    if (PAPI_query_named_event(name.c_str()) != PAPI_OK) {
                        std::cerr << "Error: PAPI event '" << name << "' not available\n";
                        return {};
                    }
                }
            }
        }
#endif
        return settings;
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

void write_json(pq_type const& pq, Settings const& settings, std::vector<ThreadData> const& data, std::ostream& out) {
    out << '{';
    out << std::quoted("pq") << ':';
    pq.write_json(out);
    out << ',' << std::quoted("settings") << ':';
    out << '{' << std::quoted("num_threads") << ':' << settings.num_threads;
    out << ',' << std::quoted("prefill_per_thread") << ':' << settings.prefill_per_thread;
    out << ',' << std::quoted("iterations_per_thread") << ':' << settings.iterations_per_thread;
    out << ',' << std::quoted("prefill_min") << ':' << settings.min_prefill;
    out << ',' << std::quoted("prefill_max") << ':' << settings.max_prefill;
    auto mode_name = [](Settings::Mode mode) {
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
    out << ',' << std::quoted("mode") << ':' << std::quoted(mode_name(settings.mode));
    if (settings.mode == Settings::Mode::Update) {
        out << ',' << std::quoted("update_min") << ':' << settings.min_update;
        out << ',' << std::quoted("update_max") << ':' << settings.max_update;
    } else if (settings.mode == Settings::Mode::Random) {
        out << ',' << std::quoted("random_min") << ':' << settings.min_random;
        out << ',' << std::quoted("random_max") << ':' << settings.max_random;
    }
    out << ',' << std::quoted("batch_size") << ':' << settings.batch_size;
    out << ',' << std::quoted("seed") << ':' << settings.seed;
#ifdef WITH_PAPI
    out << ',' << std::quoted("papi_events") << ':' << '[';
    for (std::size_t i = 0; i < settings.papi_events.size(); ++i) {
        out << std::quoted(settings.papi_events[i]);
        if (i != settings.papi_events.size() - 1) {
            out << ',';
        }
    }
    out << ']';
#endif
    out << '}';
    out << ',' << std::quoted("results") << ':';
    out << '{';
    auto last_end = std::max_element(data.begin(), data.end(), [](auto const& lhs, auto const& rhs) {
                        return lhs.end_time < rhs.end_time;
                    })->end_time;
    auto first_start = std::min_element(data.begin(), data.end(), [](auto const& lhs, auto const& rhs) {
                           return lhs.start_time < rhs.start_time;
                       })->start_time;
    out << std::quoted("total_time") << ':' << (last_end - first_start).count();
    out << ',' << std::quoted("thread_stats") << ':' << '[';
    for (auto it = data.begin(); it != data.end(); ++it) {
        out << '{';
        out << std::quoted("time") << ':' << (it->end_time - it->start_time).count();
        out << ',' << std::quoted("pushes") << ':' << it->push_count;
        out << ',' << std::quoted("pops") << ':' << it->pop_count;
        out << ',' << std::quoted("failed_pops") << ':' << it->failed_pop_count;
#ifdef WITH_PAPI
        out << ',' << std::quoted("papi_events") << ':';
        if (it->papi_success) {
            out << '[';
            for (std::size_t i = 0; i < settings.papi_events.size(); ++i) {
                out << '{';
                out << std::quoted("name") << ':' << std::quoted(settings.papi_events[i]);
                out << ',' << std::quoted("count") << ':' << it->papi_counter[i];
                out << '}';
                if (i != settings.papi_events.size() - 1) {
                    out << ',';
                }
            }
            out << ']';
        } else {
            out << "null";
        }
#endif
        out << '}';
        if (it != std::prev(data.end())) {
            out << ',';
        }
    }
    out << ']';
    out << '}';
    out << '}';
    out << '\n';
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
        std::clog << "Running...\n";
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
    context.synchronize();
    if (context.id() == 0) {
        std::clog << "Finished\n";
    }
    return std::move(context.data());
}

#ifdef LOG_OPERATIONS
void write_log(std::vector<ThreadData> const& thread_data, std::ostream& out) {
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
            out << i << ' ' << pushes[i].element.first << '\n';
            ++i;
        }
        out << '-' << ref_index << '\n';
    }
    for (; i < pushes.size(); ++i) {
        out << i << ' ' << pushes[i].element.first << '\n';
    }
}
#endif

bool run_benchmark(cxxopts::ParseResult const& parse_result) {
    auto settings = Settings::from_cmd(parse_result);
    if (!settings) {
        return false;
    }

    std::vector<ThreadData> thread_data(static_cast<std::size_t>(settings->num_threads));
    std::unique_ptr<pq_type> pq;
    auto dispatch = [&](auto& task) {
        pq = std::make_unique<pq_type>(settings->num_threads,
                                       static_cast<std::size_t>(settings->prefill_per_thread * settings->num_threads) +
                                           task.required_capacity(*settings),
                                       parse_result);
        auto dispatcher = thread_coordination::Dispatcher(
            affinity::NUMA{cores_per_numa_node, num_numa_nodes}, settings->num_threads, [&](auto t_ctx) {
                thread_data[static_cast<std::size_t>(t_ctx.id())] =
                    benchmark_thread(std::move(t_ctx), *settings, *pq, task);
            });
        dispatcher.wait();
    };

    switch (settings->mode) {
        case Settings::Mode::Update: {
            auto task = Update(*settings);
            dispatch(task);
            break;
        }
        case Settings::Mode::Random: {
            auto task = Random(*settings);
            dispatch(task);
            break;
        }
        case Settings::Mode::PopAll: {
            auto task = PopAll();
            dispatch(task);
            break;
        }
        case Settings::Mode::PushPop: {
            auto task = PushPop(*settings);
            dispatch(task);
            break;
        }
    }

    write_json(*pq, *settings, thread_data, std::cout);
#ifdef LOG_OPERATIONS
    write_log(thread_data, settings->log_out);
#endif
    return true;
}

int main(int argc, char* argv[]) {
    write_build_info(std::clog);
    std::clog << "\nCommand line:";
    for (int i = 0; i < argc; ++i) {
        std::clog << ' ' << argv[i];
    }
    std::clog << '\n';

    cxxopts::Options cmd(argv[0]);
    cmd.add_options()("h,help", "Print this help", cxxopts::value<bool>());
    Settings::add_cmd_options(cmd);
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

    bool success = run_benchmark(args);
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
