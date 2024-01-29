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
#include <new>
#include <numeric>
#include <random>
#include <type_traits>
#include <utility>
#include <vector>

using key_type = unsigned long;
using value_type = unsigned long;

using pq_type = PQ<true, key_type, value_type>;
using handle_type = pq_type::handle_type;

struct Settings {
    int num_threads = 4;
    long long prefill_per_thread = 1 << 20;
    long long iterations_per_thread = 1 << 22;
    key_type min_prefill = 1;
    key_type max_prefill = 1 << 20;
    long min_update = 0;
    long max_update = 1 << 20;
    long long batch_size = 1 << 12;
    int seed = 1;
#ifdef LOG_OPERATIONS
    std::filesystem::path log_file = "operation_log.txt";
#endif
#ifdef WITH_PAPI
    std::vector<std::string> papi_events;
#endif
    pq_type::settings_type pq_settings;
};

void register_cmd_options(Settings& settings, cxxopts::Options& cmd) {
    cmd.add_options()
        // clang-format off
            ("j,threads", "Number of threads", cxxopts::value<int>(settings.num_threads), "NUMBER")
            ("p,prefill", "Prefill per thread", cxxopts::value<long long>(settings.prefill_per_thread), "NUMBER")
            ("n,iterations", "Number of iterations per thread", cxxopts::value<long long>(settings.iterations_per_thread), "NUMBER")
            ("min-prefill", "Min prefill key", cxxopts::value<key_type>(settings.min_prefill), "NUMBER")
            ("max-prefill", "Max prefill key", cxxopts::value<key_type>(settings.max_prefill), "NUMBER")
            ("min-update", "Min update", cxxopts::value<long>(settings.min_update), "NUMBER")
            ("max-update", "Max update", cxxopts::value<long>(settings.max_update), "NUMBER")
            ("batch-size", "Batch size", cxxopts::value<long long>(settings.batch_size), "NUMBER")
            ("s,seed", "Initial seed", cxxopts::value<int>(settings.seed), "NUMBER")
#ifdef LOG_OPERATIONS
            ("l,log-file", "File to write the operation log to", cxxopts::value<std::filesystem::path>(settings.log_file), "PATH")
#endif
#ifdef WITH_PAPI
            ("r,count-event", "Papi event to count", cxxopts::value<std::vector<std::string>>(settings.papi_events), "STRING")
#endif
        // clang-format on
        ;
    settings.pq_settings.register_cmd_options(cmd);
}

bool validate_settings(Settings const& settings) {
    if (settings.num_threads <= 0) {
        std::cerr << "Error: Number of threads must be greater than 0\n";
        return false;
    }
    if (settings.prefill_per_thread < 0) {
        std::cerr << "Error: Prefill must be nonnegative\n";
        return false;
    }
    if (settings.iterations_per_thread < 0) {
        std::cerr << "Error: Iterations must be nonnegative\n";
        return false;
    }
    if (settings.min_prefill <= 0) {
        std::cerr << "Error: Prefill keys must be greater than 0\n";
        return false;
    }
    if (settings.max_prefill < settings.min_prefill) {
        std::cerr << "Error: Invalid prefill range\n";
        return false;
    }
    if (settings.min_update < 0) {
        std::cerr << "Error: Min update must be nonnegative\n";
        return false;
    }
    if (settings.max_update < settings.min_update) {
        std::cerr << "Error: Invalid update range\n";
        return false;
    }
    if (settings.batch_size <= 0) {
        std::cerr << "Error: batch size must be greater than 0\n";
        return false;
    }
    if (settings.seed <= 0) {
        std::cerr << "Error: Seed must be greater than 0\n";
        return false;
    }
#ifdef LOG_OPERATIONS
    if (settings.log_file.empty()) {
        std::cerr << "Error: Log file name must not be empty\n";
        return false;
    }
    auto out = std::ofstream(settings.log_file);
    if (out.fail()) {
        std::cerr << "Error: Could not open file " << settings.log_file << " for writing\n";
        return false;
    }
    out.close();
#endif
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
    return settings.pq_settings.validate();
}

void write_settings_human_readable(Settings const& settings, std::ostream& out) {
    out << "Threads: " << settings.num_threads << '\n';
    out << "Prefill per thread: " << settings.prefill_per_thread << '\n';
    out << "Iterations per thread: " << settings.iterations_per_thread << '\n';
    out << "Prefill range: [" << settings.min_prefill << ", " << settings.max_prefill << "]\n";
    out << "Update range: [" << settings.min_update << ", " << settings.max_update << "]\n";
    out << "Batch size: " << settings.batch_size << '\n';
    out << "Seed: " << settings.seed << '\n';
#ifdef LOG_OPERATIONS
    out << "Log file: " << settings.log_file << '\n';
#endif
#ifdef WITH_PAPI
    out << "PAPI events: [";
    for (std::size_t i = 0; i < settings.papi_events.size(); ++i) {
        out << settings.papi_events[i];
        if (i != settings.papi_events.size() - 1) {
            out << ", ";
        }
    }
    out << ']' << '\n';
#endif
    settings.pq_settings.write_human_readable(out);
}

void write_settings_json(Settings const& settings, std::ostream& out) {
    out << '{';
    out << std::quoted("num-threads") << ':' << settings.num_threads << ',';
    out << std::quoted("prefill-per-thread") << ':' << settings.prefill_per_thread << ',';
    out << std::quoted("iterations-per-thread") << ':' << settings.iterations_per_thread << ',';
    out << std::quoted("prefill-min") << ':' << settings.min_prefill << ',';
    out << std::quoted("prefill-max") << ':' << settings.max_prefill << ',';
    out << std::quoted("update-min") << ':' << settings.min_update << ',';
    out << std::quoted("update-max") << ':' << settings.max_update << ',';
    out << std::quoted("batch-size") << ':' << settings.batch_size << ',';
    out << std::quoted("seed") << ':' << settings.seed << ',';
#ifdef WITH_PAPI
    out << std::quoted("papi-events") << ':';
    out << '[';
    for (std::size_t i = 0; i < settings.papi_events.size(); ++i) {
        out << std::quoted(settings.papi_events[i]);
        if (i != settings.papi_events.size() - 1) {
            out << ',';
        }
    }
    out << ']' << ',';
#endif
    out << std::quoted("pq") << ':';
    settings.pq_settings.write_json(out);
    out << '}';
}

struct ThreadData {
    std::chrono::high_resolution_clock::time_point start_time{};
    std::chrono::high_resolution_clock::time_point end_time{};
    long long iter_count = 0;
    long long failed_pop_count = 0;
#ifdef WITH_PAPI
    std::vector<long long> papi_event_counter{};
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

void write_thread_data_json(ThreadData const& data, std::ostream& out) {
    out << '{';
    out << std::quoted("time-ns") << ':' << std::chrono::nanoseconds{data.end_time - data.start_time}.count() << ',';
    out << std::quoted("iterations") << ':' << data.iter_count << ',';
    out << std::quoted("failed-pops") << ':' << data.failed_pop_count;
#ifdef WITH_PAPI
    out << ',';
    out << std::quoted("papi-event-counter") << ':';
    out << '[';
    for (std::size_t i = 0; i < data.papi_event_counter.size(); ++i) {
        out << data.papi_event_counter[i];
        if (i != data.papi_event_counter.size() - 1) {
            out << ',';
        }
    }
    out << ']';
#endif
    out << '}';
}

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

    ThreadData& data() noexcept {
        return data_;
    }
    ThreadData const& data() const noexcept {
        return data_;
    }
};

class Task {
    std::vector<long long> data_;
    std::atomic_llong counter_{0};

   public:
    explicit Task(Settings const& settings)
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
            context.data().iter_count += to - from;
        }
    }
};

void write_result_json(Settings const& settings, std::vector<ThreadData> const& data, std::ostream& out) {
    out << '{';
    out << std::quoted("settings") << ':';
    write_settings_json(settings, out);
    out << ',';
    auto last_end = std::max_element(data.begin(), data.end(), [](auto const& lhs, auto const& rhs) {
                        return lhs.end_time < rhs.end_time;
                    })->end_time;
    auto first_start = std::min_element(data.begin(), data.end(), [](auto const& lhs, auto const& rhs) {
                           return lhs.start_time < rhs.start_time;
                       })->start_time;
    out << std::quoted("time-ns") << ':' << std::chrono::nanoseconds{last_end - first_start}.count() << ',';
    out << std::quoted("per-thread") << ':';
    out << '[';
    for (auto it = data.begin(); it != data.end(); ++it) {
        write_thread_data_json(*it, out);
        if (it != std::prev(data.end())) {
            out << ',';
        }
    }
    out << ']';
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

ThreadData benchmark_thread(thread_coordination::Context thread_context, Settings const& settings, pq_type& pq,
                            Task& task) {
    Context context{std::move(thread_context), pq.get_handle()};
#ifdef WITH_PAPI
    int event_set = PAPI_NULL;
    if (!settings.papi_events.empty()) {
        context.data().papi_event_counter.resize(settings.papi_events.size());
        try {
            event_set = prepare_papi(settings);
        } catch (std::exception const& e) {
            context.write(std::cerr) << e.what() << '\n';
        }
    }
#endif
#ifdef LOG_OPERATIONS
    context.data().pushes.reserve(
        static_cast<std::size_t>(settings.prefill_per_thread + 2 * settings.iterations_per_thread));
    context.data().pops.reserve(static_cast<std::size_t>(2 * settings.iterations_per_thread));
#endif

    std::vector<key_type> prefill(static_cast<std::size_t>(settings.prefill_per_thread));

    if (context.id() == 0) {
        std::clog << "Preparing...\n";
    }
    std::seed_seq seed{settings.seed, context.id()};
    std::default_random_engine rng(seed);
    context.synchronize();
    std::generate(prefill.begin(), prefill.end(), [&rng, min = settings.min_prefill, max = settings.max_prefill]() {
        return std::uniform_int_distribution<key_type>(min, max)(rng);
    });
    context.synchronize();
    task.prepare(settings, context, rng);
    context.synchronize();
    if (context.id() == 0) {
        std::clog << "Prefilling...\n";
    }
    context.synchronize();
    for (auto i = 0LL; i < settings.prefill_per_thread; ++i) {
        context.push({prefill[static_cast<std::size_t>(i)],
                      static_cast<value_type>(context.id() * settings.prefill_per_thread + i)});
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
        if (int ret = PAPI_stop(event_set, context.data().papi_event_counter.data()); ret != PAPI_OK) {
            context.write(std::cerr) << "Failed to stop performance counters\n";
        }
    }
#endif
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
        push_index[pushes[i].element.second] = i;
    }
    std::sort(pops.begin(), pops.end(), [](auto const& lhs, auto const& rhs) { return lhs.tick < rhs.tick; });
    out << pushes.size() << ' ' << pops.size() << '\n';
    std::size_t i = 0;
    for (auto const& pop : pops) {
        while ((i != pushes.size() && pushes[i].tick < pop.tick)) {
            out << '+' << pushes[i].element.first << '\n';
            ++i;
        }
        out << '-' << push_index[static_cast<std::size_t>(pop.val)] << '\n';
    }
    for (; i < pushes.size(); ++i) {
        out << i << ' ' << pushes[i].element.first << '\n';
    }
}
#endif

void run_benchmark(Settings const& settings) {
    Task task(settings);
    std::vector<ThreadData> thread_data(static_cast<std::size_t>(settings.num_threads));
    auto pq =
        pq_type(settings.num_threads, static_cast<std::size_t>(settings.prefill_per_thread * settings.num_threads),
                settings.pq_settings);

    auto dispatcher = thread_coordination::Dispatcher(settings.num_threads, [&](auto ctx) {
        auto t_id = static_cast<std::size_t>(ctx.id());
        thread_data[t_id] = benchmark_thread(std::move(ctx), settings, pq, task);
    });
    dispatcher.wait();

#ifdef LOG_OPERATIONS
    std::clog << "Writing logs...\n";
    std::ofstream log_out(settings.log_file);  // assumed to be valid
    write_log(thread_data, log_out);
    log_out.close();
#endif
    std::clog << "Done\n";
    write_result_json(settings, thread_data, std::cout);
}

int main(int argc, char* argv[]) {
    write_build_info(std::clog);
    std::clog << '\n';

    std::clog << "= Priority queue =\n";
    pq_type::write_human_readable(std::clog);
    std::clog << '\n';

    std::clog << "= Command line =\n";
    for (int i = 0; i < argc; ++i) {
        std::clog << argv[i];
        if (i != argc - 1) {
            std::clog << ' ';
        }
    }
    std::clog << '\n' << '\n';

    cxxopts::Options cmd(argv[0]);
    cmd.add_options()("h,help", "Print this help", cxxopts::value<bool>());
    Settings settings{};
    register_cmd_options(settings, cmd);

    try {
        auto args = cmd.parse(argc, argv);
        if (args.count("help") > 0) {
            std::clog << cmd.help() << '\n';
            return EXIT_SUCCESS;
        }
    } catch (std::exception const& e) {
        std::cerr << "Error parsing command line: " << e.what() << '\n';
        std::cerr << "Use --help for usage information" << '\n';
        return EXIT_FAILURE;
    }

    std::clog << "= Settings =\n";
    write_settings_human_readable(settings, std::clog);
    std::clog << '\n';

    if (!validate_settings(settings)) {
        return EXIT_FAILURE;
    }

    std::clog << "= Running benchmark =\n";
    run_benchmark(settings);
    return EXIT_SUCCESS;
}
