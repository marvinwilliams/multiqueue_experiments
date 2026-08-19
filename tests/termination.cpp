#include "util/build_info.hpp"
#include "util/selector.hpp"
#include "util/termination_detection.hpp"
#include "util/thread_coordination.hpp"

#include <cxxopts.hpp>

#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <random>
#include <utility>

using key_type = unsigned long;
using value_type = unsigned long;

using pq_type = PQ<true, key_type, value_type>;
using handle_type = pq_type::handle_type;

struct Settings {
    int num_threads = 4;
    long long prefill_per_thread = 1 << 20;
    key_type min_prefill = 1;
    key_type max_prefill = 1 << 20;
    int push_probability_percent = 0;
    int seed = 1;
#ifdef LOG_OPERATIONS
    std::filesystem::path log_file = "operation_log.txt";
#endif
    pq_type::settings_type pq_settings;
};

Settings settings{};

void register_cmd_options(cxxopts::Options& cmd) {
    cmd.add_options()
        // clang-format off
            ("j,threads", "Number of threads", cxxopts::value<int>(settings.num_threads), "NUMBER")
            ("p,prefill", "Prefill per thread", cxxopts::value<long long>(settings.prefill_per_thread), "NUMBER")
            ("min-prefill", "Min prefill key", cxxopts::value<key_type>(settings.min_prefill), "NUMBER")
            ("max-prefill", "Max prefill key", cxxopts::value<key_type>(settings.max_prefill), "NUMBER")
            ("r,push_probability_percent", "Probability to push an element in percent", cxxopts::value<int>(settings.push_probability_percent), "NUMBER")
            ("s,seed", "Initial seed", cxxopts::value<int>(settings.seed), "NUMBER")
#ifdef LOG_OPERATIONS
            ("l,log-file", "File to write the operation log to", cxxopts::value<std::filesystem::path>(settings.log_file), "PATH")
#endif
        // clang-format on
        ;
    settings.pq_settings.register_cmd_options(cmd);
}

bool validate_settings() {
    if (settings.num_threads <= 0) {
        std::cerr << "Error: Number of threads must be greater than 0\n";
        return false;
    }
    if (settings.prefill_per_thread < 0) {
        std::cerr << "Error: Prefill must be nonnegative\n";
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
    if (settings.push_probability_percent < 0 || settings.push_probability_percent > 100) {
        std::cerr << "Error: Push rate must be between 0 and 100\n";
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
    return settings.pq_settings.validate();
}

void write_settings_human_readable(std::ostream& out) {
    out << "Threads: " << settings.num_threads << '\n';
    out << "Prefill per thread: " << settings.prefill_per_thread << '\n';
    out << "Prefill range: [" << settings.min_prefill << ", " << settings.max_prefill << "]\n";
    out << "Push rate: " << settings.push_probability_percent << "%\n";
    out << "Seed: " << settings.seed << '\n';
#ifdef LOG_OPERATIONS
    out << "Log file: " << settings.log_file << '\n';
#endif
    settings.pq_settings.write_human_readable(out);
}

void write_settings_json(std::ostream& out) {
    out << '{';
    out << std::quoted("num_threads") << ':' << settings.num_threads << ',';
    out << std::quoted("prefill_per_thread") << ':' << settings.prefill_per_thread << ',';
    out << std::quoted("prefill_min") << ':' << settings.min_prefill << ',';
    out << std::quoted("prefill_max") << ':' << settings.max_prefill << ',';
    out << std::quoted("push_probability_percent") << ':' << settings.push_probability_percent << ',';
    out << std::quoted("seed") << ':' << settings.seed << ',';
    out << std::quoted("pq") << ':';
    settings.pq_settings.write_json(out);
    out << '}';
}

struct ThreadData {
    long long pop_count = 0;
    long long failed_pop_count = 0;
    std::default_random_engine rng;
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
    out << std::quoted("pops") << ':' << data.pop_count << ',';
    out << std::quoted("failed_pops") << ':' << data.failed_pop_count;
    out << '}';
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

struct SharedData {
    std::vector<long long> updates;
    termination_detection::TerminationDetection termination_detection;
    std::vector<char> found;
    std::atomic_bool done{false};
    std::atomic_llong missing_nodes{0};
    std::chrono::high_resolution_clock::time_point start_time;
    std::chrono::high_resolution_clock::time_point end_time;
    std::vector<ThreadData> thread_data;
};

SharedData shared_data;

void write_result_json(std::ostream& out) {
    out << '{';
    out << std::quoted("settings") << ':';
    write_settings_json(out);
    out << ',';
    out << std::quoted("results") << ':';
    out << '{';
    out << std::quoted("time_ns") << ':'
        << std::chrono::nanoseconds{shared_data.end_time - shared_data.start_time}.count() << ',';
    out << std::quoted("thread_data") << ':';
    out << '[';
    for (auto it = shared_data.thread_data.begin(); it != shared_data.thread_data.end(); ++it) {
        write_thread_data_json(*it, out);
        if (it != std::prev(shared_data.thread_data.end())) {
            out << ',';
        }
    }
    out << ']';
    out << '}';
    out << '}' << '\n';
}

class Context : public thread_coordination::Context {
    handle_type handle_;
    ThreadData thread_data_;

   public:
    explicit Context(thread_coordination::Context ctx, handle_type handle)
        : thread_coordination::Context{std::move(ctx)}, handle_{std::move(handle)} {
    }
#ifdef LOG_OPERATIONS
    void push(std::pair<key_type, value_type> const& e) {
        handle_.push(e);
        auto tick = std::chrono::high_resolution_clock::now();
        thread_data_.pushes.push_back({tick, e});
    }

    auto try_pop() {
        auto tick = std::chrono::high_resolution_clock::now();
        auto retval = handle_.try_pop();
        if (retval) {
            thread_data_.pops.push_back({tick, retval->second});
        }
        return retval;
    }
#else
    void push(std::pair<key_type, value_type> const& e) {
        handle_.push(e);
    }

    auto try_pop() {
        return handle_.try_pop();
    }
#endif

    ThreadData& thread_data() noexcept {
        return thread_data_;
    }
    [[nodiscard]] ThreadData const& thread_data() const noexcept {
        return thread_data_;
    }
};

[[gnu::noinline]] void work_loop(Context& context) {
    auto push_dist = std::bernoulli_distribution{settings.push_probability_percent / 100.0};
    while (true) {
        context.synchronize();
        if (context.id() == 0) {
            shared_data.termination_detection.reset();
            shared_data.done.store(true, std::memory_order_relaxed);
        }
        auto pop_count = 0;
        auto push_count = 0;
        context.synchronize();
        while (shared_data.termination_detection.repeat([&]() {
            auto node = context.try_pop();
            if (!node) {
                ++context.thread_data().failed_pop_count;
                return false;
            }
            ++shared_data.found[static_cast<std::size_t>(node->first)];
            ++pop_count;
            return true;
        })) {
            if (settings.push_probability_percent > 0 && push_dist(context.thread_data().rng)) {
                auto dist = std::uniform_int_distribution<value_type>(settings.min_prefill, settings.max_prefill);
                context.push({dist(context.thread_data().rng), dist(context.thread_data().rng)});
                ++push_count;
            }
        }
        shared_data.missing_nodes.fetch_add(push_count - pop_count, std::memory_order_relaxed);
        context.thread_data().pop_count += pop_count;
        context.synchronize();
        if (context.id() == 0) {
            std::clog << "Missing nodes: " << shared_data.missing_nodes.load(std::memory_order_relaxed) << '\n';
            auto unpopped =
                std::count_if(shared_data.found.begin(), shared_data.found.end(), [](auto e) { return e == 0; });
            auto popped_twice =
                std::count_if(shared_data.found.begin(), shared_data.found.end(), [](auto e) { return e > 1; });
            std::clog << "Unpopped nodes: " << unpopped << '\n';
            std::clog << "Popped twice: " << popped_twice << '\n';
        }
        context.synchronize();
        auto node = context.try_pop();
        if (node) {
            ++shared_data.found[static_cast<std::size_t>(node->first)];
            ++context.thread_data().pop_count;
            shared_data.done.store(false, std::memory_order_relaxed);
            shared_data.missing_nodes.fetch_sub(1, std::memory_order_relaxed);
        }
        context.synchronize();
        if (context.id() == 0) {
            std::clog << "Missing nodes after synchronization: "
                      << shared_data.missing_nodes.load(std::memory_order_relaxed) << '\n';
            auto unpopped =
                std::count_if(shared_data.found.begin(), shared_data.found.end(), [](auto e) { return e == 0; });
            auto popped_twice =
                std::count_if(shared_data.found.begin(), shared_data.found.end(), [](auto e) { return e > 1; });
            std::clog << "Unpopped nodes: " << unpopped << '\n';
            std::clog << "Popped twice: " << popped_twice << '\n';
        }
        context.synchronize();
        if (shared_data.missing_nodes.load(std::memory_order_relaxed) == 0) {
            break;
        }
    }
}

void benchmark_thread(Context context) {
#ifdef LOG_OPERATIONS
    context.thread_data().pushes.reserve(static_cast<std::size_t>(2 * settings.prefill_per_thread));
    context.thread_data().pops.reserve(static_cast<std::size_t>(2 * settings.prefill_per_thread));
#endif

    std::vector<key_type> prefill(static_cast<std::size_t>(settings.prefill_per_thread));

    if (context.id() == 0) {
        std::clog << "Preparing...\n";
    }
    std::seed_seq seed{settings.seed, context.id()};
    context.thread_data().rng.seed(seed);
    context.synchronize();
    std::generate(prefill.begin(), prefill.end(), [&, min = settings.min_prefill, max = settings.max_prefill]() {
        return std::uniform_int_distribution<key_type>(min, max)(context.thread_data().rng);
    });
    std::iota(prefill.begin(), prefill.end(), static_cast<key_type>(context.id() * settings.prefill_per_thread));
    std::shuffle(prefill.begin(), prefill.end(), context.thread_data().rng);
    context.synchronize();
    if (context.id() == 0) {
        std::clog << "Prefilling...\n";
    }
    context.synchronize();
    for (auto i = 0LL; i < settings.prefill_per_thread; ++i) {
        context.push({prefill[static_cast<std::size_t>(i)],
                      static_cast<value_type>(context.id() * settings.prefill_per_thread + i)});
        /* context.push({prefill[static_cast<std::size_t>(i)], 1}); */
    }
    context.synchronize();
    if (context.id() == 0) {
        std::clog << "Working...\n";
    }
    if (context.id() == 0) {
        shared_data.start_time = std::chrono::high_resolution_clock::now();
    }
    context.synchronize();
    work_loop(context);
    context.synchronize();
    if (context.id() == 0) {
        shared_data.end_time = std::chrono::high_resolution_clock::now();
    }
    shared_data.thread_data[static_cast<std::size_t>(context.id())] = std::move(context.thread_data());
}

void run_benchmark() {
    shared_data.thread_data.resize(static_cast<std::size_t>(settings.num_threads));
    shared_data.termination_detection.reset(settings.num_threads);
    shared_data.missing_nodes.store(settings.num_threads * settings.prefill_per_thread, std::memory_order_relaxed);
    shared_data.found.resize(settings.num_threads * settings.prefill_per_thread, 0);
    auto pq =
        pq_type(settings.num_threads, static_cast<std::size_t>(settings.prefill_per_thread * settings.num_threads),
                settings.pq_settings);

    auto dispatcher =
        thread_coordination::Dispatcher(thread_coordination::affinity::None{}, settings.num_threads,
                                        [&](auto ctx) { benchmark_thread(Context(std::move(ctx), pq.get_handle())); });
    dispatcher.wait();

#ifdef LOG_OPERATIONS
    std::clog << "Writing logs...\n";
    std::ofstream log_out(settings.log_file);  // assumed to be valid
    write_log(shared_data.thread_data, log_out);
    log_out.close();
#endif
    std::clog << "Done\n";
    std::clog << '\n';
    std::clog << "= Results =\n";
    std::clog << "Time (s): " << std::fixed << std::setprecision(3)
              << std::chrono::duration<double>(shared_data.end_time - shared_data.start_time).count() << '\n';
    write_result_json(std::cout);
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
    register_cmd_options(cmd);

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
    write_settings_human_readable(std::clog);
    std::clog << '\n';

    if (!validate_settings()) {
        return EXIT_FAILURE;
    }

    std::clog << "= Running benchmark =\n";
    run_benchmark();
    return EXIT_SUCCESS;
}
