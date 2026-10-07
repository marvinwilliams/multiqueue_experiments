#pragma once

#include <pthread.h>
#include <sched.h>

#ifdef __SSE2__
#include <emmintrin.h>
#endif

#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cctype>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <numeric>
#include <ostream>
#include <string>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

namespace thread_coordination {

inline void cpu_relax() noexcept {
#if defined(__SSE2__)
    _mm_pause();
#elif defined(__aarch64__)
    asm volatile("isb" ::: "memory");
#endif
}

namespace affinity {

enum Policy : int {
    None = 0,
    ThreadId = 1,
    Same = 2,
    CloseCaches = 3,
    FarCaches = 4,
    CloseL3FarL1 = 5,
    FarL1CloseL3 = 6,
    CoresFirst = 7,
};

inline constexpr int max_id = CoresFirst;
inline constexpr int default_id = CoresFirst;

inline char const* name(int id) noexcept {
    switch (id) {
        case None:
            return "None";
        case ThreadId:
            return "Thread Id";
        case Same:
            return "Same";
        case CloseCaches:
            return "Close caches";
        case FarCaches:
            return "Far caches";
        case CloseL3FarL1:
            return "Close L3 Far L1";
        case FarL1CloseL3:
            return "Far L1 Close L3";
        case CoresFirst:
            return "Cores first (L3, NUMA node, then SMT)";
        default:
            return "Unknown";
    }
}

inline std::vector<std::size_t> available_cpus() {
    cpu_set_t set;
    CPU_ZERO(&set);
    std::vector<std::size_t> cpus;
    if (sched_getaffinity(0, sizeof(set), &set) != 0) {
        cpus.resize(std::thread::hardware_concurrency());
        std::iota(cpus.begin(), cpus.end(), 0);
        return cpus;
    }
    for (std::size_t i = 0; i < CPU_SETSIZE; ++i) {
        if (CPU_ISSET(i, &set)) {
            cpus.push_back(i);
        }
    }
    return cpus;
}

inline bool needs_distinct_cpus(int id) noexcept {
    return id != None && id != Same;
}

namespace detail {

struct CpuInfo {
    std::size_t cpu;
    // [0]: rank of the CPU within its core (SMT sibling)
    // [1]: rank of the core within its L3
    // [2]: rank of the L3 within its NUMA node
    // [3]: rank of the NUMA node within its package
    // [4]: id of the package
    std::array<std::size_t, 5> level;
};

inline std::size_t numa_node(std::size_t cpu) {
    std::error_code ec;
    auto dir = std::filesystem::directory_iterator("/sys/devices/system/cpu/cpu" + std::to_string(cpu), ec);
    for (; !ec && dir != std::filesystem::directory_iterator{}; dir.increment(ec)) {
        auto name = dir->path().filename().string();
        if (name.size() > 4 && name.compare(0, 4, "node") == 0 &&
            std::all_of(name.begin() + 4, name.end(), [](unsigned char c) { return std::isdigit(c) != 0; })) {
            return std::strtoul(name.c_str() + 4, nullptr, 10);
        }
    }
    return 0;
}

inline std::size_t first_cpu_of(std::string const& list_file, std::size_t fallback) {
    std::string list;
    std::ifstream(list_file) >> list;
    auto end = list.find_first_not_of("0123456789");
    return end == 0 || list.empty() ? fallback : std::strtoul(list.c_str(), nullptr, 10);
}

inline std::size_t l3_group(std::size_t cpu) {
    auto base = "/sys/devices/system/cpu/cpu" + std::to_string(cpu) + "/cache/";
    for (int i = 0;; ++i) {
        auto index = base + "index" + std::to_string(i);
        int level = 0;
        if (!(std::ifstream(index + "/level") >> level)) {
            return 0;
        }
        if (level == 3) {
            return first_cpu_of(index + "/shared_cpu_list", 0);
        }
    }
}

inline std::vector<CpuInfo> read_topology(std::vector<std::size_t> const& cpus) {
    std::vector<CpuInfo> hierarchy(cpus.size());
    for (std::size_t i = 0; i < cpus.size(); ++i) {
        auto cpu = cpus[i];
        auto topology = "/sys/devices/system/cpu/cpu" + std::to_string(cpu) + "/topology/";
        auto& h = hierarchy[i].level;
        hierarchy[i].cpu = cpu;
        h[0] = cpu;
        h[1] = first_cpu_of(topology + "thread_siblings_list", cpu);
        h[2] = l3_group(cpu);
        h[3] = numa_node(cpu);
        h[4] = first_cpu_of(topology + "physical_package_id", 0);
    }
    return hierarchy;
}

inline void rank_levels(std::vector<CpuInfo>& hierarchy) {
    for (std::size_t l = 0; l + 1 < std::tuple_size_v<decltype(CpuInfo::level)>; ++l) {
        std::vector<std::vector<std::size_t>> lookup;
        for (auto& info : hierarchy) {
            auto& h = info.level;
            if (h[l + 1] >= lookup.size()) {
                lookup.resize(h[l + 1] + 1);
            }
            auto& group = lookup[h[l + 1]];
            auto it = std::find(group.begin(), group.end(), h[l]);
            if (it == group.end()) {
                group.push_back(h[l]);
                it = group.end() - 1;
            }
            h[l] = static_cast<std::size_t>(std::distance(group.begin(), it));
        }
    }
}

inline std::vector<std::size_t> order_by(std::vector<CpuInfo> hierarchy, std::array<std::size_t, 5> levels) {
    rank_levels(hierarchy);
    auto key = [&levels](CpuInfo const& c) {
        return std::tie(c.level[levels[0]], c.level[levels[1]], c.level[levels[2]], c.level[levels[3]],
                        c.level[levels[4]]);
    };
    std::sort(hierarchy.begin(), hierarchy.end(),
              [&key](CpuInfo const& a, CpuInfo const& b) { return key(a) < key(b); });
    std::vector<std::size_t> order(hierarchy.size());
    std::transform(hierarchy.begin(), hierarchy.end(), order.begin(), [](CpuInfo const& c) { return c.cpu; });
    return order;
}

}  // namespace detail

using detail::numa_node;

inline std::vector<std::size_t> cpu_assignment(int id, int num_threads) {
    auto const n = static_cast<std::size_t>(num_threads);
    if (id == None) {
        return {};
    }
    auto cpus = available_cpus();
    assert(!cpus.empty());
    if (id == Same) {
        return std::vector<std::size_t>(n, cpus[0]);
    }
    assert(n <= cpus.size());
    std::vector<std::size_t> order;
    switch (id) {
        case ThreadId:
            order = std::move(cpus);
            break;
        case CloseCaches:
            order = detail::order_by(detail::read_topology(cpus), {4, 3, 2, 1, 0});
            break;
        case FarCaches:
            order = detail::order_by(detail::read_topology(cpus), {4, 0, 1, 2, 3});
            break;
        case CloseL3FarL1:
            order = detail::order_by(detail::read_topology(cpus), {4, 3, 2, 0, 1});
            break;
        case FarL1CloseL3:
            order = detail::order_by(detail::read_topology(cpus), {4, 0, 3, 2, 1});
            break;
        case CoresFirst:
            order = detail::order_by(detail::read_topology(cpus), {0, 4, 3, 2, 1});
            break;
        default:
            assert(false && "invalid affinity");
    }
    order.resize(n);
    return order;
}

}  // namespace affinity

class Barrier {
    pthread_barrier_t barrier_{};

   public:
    explicit Barrier(int num_threads) {
        assert(num_threads > 0);
        if (int rc = pthread_barrier_init(&barrier_, nullptr, static_cast<unsigned int>(num_threads)); rc != 0) {
            std::cerr << "Error: " << "Failed to create barrier" << ": " << std::strerror(rc) << std::endl;
            std::abort();
        }
    }

    Barrier(Barrier const&) = delete;
    Barrier& operator=(Barrier const&) = delete;
    Barrier(Barrier&&) = delete;
    Barrier& operator=(Barrier&&) = delete;

    ~Barrier() {
        pthread_barrier_destroy(&barrier_);
    }

    void wait() {
        pthread_barrier_wait(&barrier_);
    }
};

class SpinBarrier {
    int num_threads_;
    alignas(64) std::atomic_int count_{0};
    alignas(64) std::atomic_uint generation_{0};

   public:
    explicit SpinBarrier(int num_threads) noexcept : num_threads_{num_threads} {
    }

    void wait() noexcept {
        auto generation = generation_.load(std::memory_order_relaxed);
        if (count_.fetch_add(1, std::memory_order_acq_rel) + 1 == num_threads_) {
            count_.store(0, std::memory_order_relaxed);
            generation_.store(generation + 1, std::memory_order_release);
            return;
        }
        while (generation_.load(std::memory_order_relaxed) == generation) {
            cpu_relax();
        }
        std::atomic_thread_fence(std::memory_order_acquire);
    }
};

class Context {
    friend class Dispatcher;

    int id_;
    Barrier* barrier_;
    SpinBarrier* spin_barrier_;
    std::mutex* write_mutex_{};

    Context(int id, Barrier& barrier, SpinBarrier& spin_barrier, std::mutex& m)
        : id_{id}, barrier_{&barrier}, spin_barrier_{&spin_barrier}, write_mutex_{&m} {
    }

    class GuardedWriter {
        friend Context;
        std::unique_lock<std::mutex> lock_;
        std::ostream* out_;

        explicit GuardedWriter(std::mutex& mutex, std::ostream& out, int id) : lock_{mutex}, out_{&out} {
            *out_ << "[Thread " << id << "] ";
        }

       public:
        template <typename T>
        std::ostream& operator<<(T&& t) {
            *out_ << std::forward<T>(t);
            return *out_;
        }
    };

   public:
    Context(Context const&) = delete;
    Context(Context&&) noexcept = default;

    Context& operator=(Context const&) = delete;
    Context& operator=(Context&&) noexcept = default;

    ~Context() = default;

    [[nodiscard]] int id() const noexcept {
        return id_;
    }

    GuardedWriter write(std::ostream& out) {
        return GuardedWriter(*write_mutex_, out, id_);
    }

    void synchronize() const {
        barrier_->wait();
    }

    void spin_synchronize() const {
        spin_barrier_->wait();
    }
};

class Dispatcher {
    std::vector<pthread_t> threads_;
    Barrier barrier_;
    SpinBarrier spin_barrier_;
    std::mutex write_mutex_;

    template <typename F>
    static void* trampoline(void* f) {
        std::unique_ptr<F>{static_cast<F*>(f)}->operator()();
        return nullptr;
    }

    template <typename F>
    void start(F f, std::size_t const* cpu) {
        pthread_attr_t attr;
        if (int rc = pthread_attr_init(&attr); rc != 0) {
            std::cerr << "Error: " << "Failed to initialize thread attributes" << ": " << std::strerror(rc)
                      << std::endl;
            std::abort();
        }
        if (cpu != nullptr) {
            cpu_set_t set;
            CPU_ZERO(&set);
            CPU_SET(*cpu, &set);
            if (int rc = pthread_attr_setaffinity_np(&attr, sizeof(set), &set); rc != 0) {
                std::cerr << "Error: " << "Failed to pin thread to CPU " + std::to_string(*cpu) << ": "
                          << std::strerror(rc) << std::endl;
                std::abort();
            }
        }
        auto owned = std::make_unique<F>(std::move(f));
        pthread_t thread{};
        if (int rc = pthread_create(&thread, &attr, &trampoline<F>, owned.get()); rc != 0) {
            std::cerr << "Error: " << "Failed to create thread" << ": " << std::strerror(rc) << std::endl;
            std::abort();
        }
        owned.release();
        pthread_attr_destroy(&attr);
        threads_.push_back(thread);
    }

   public:
    template <typename Task, typename... Args>
    Dispatcher(std::vector<std::size_t> const& cpus, int num_threads, Task task, Args... args)
        : barrier_(num_threads), spin_barrier_(num_threads) {
        assert(cpus.empty() || cpus.size() >= static_cast<std::size_t>(num_threads));
        threads_.reserve(static_cast<std::size_t>(num_threads));
        for (int i = 0; i < num_threads; ++i) {
            start([task, ctx = Context{i, barrier_, spin_barrier_, write_mutex_},
                   args...]() mutable { task(std::move(ctx), args...); },
                  cpus.empty() ? nullptr : &cpus[static_cast<std::size_t>(i)]);
        }
    }

    Dispatcher(Dispatcher const&) = delete;
    Dispatcher& operator=(Dispatcher const&) = delete;

    void wait() {
        for (auto thread : threads_) {
            pthread_join(thread, nullptr);
        }
        threads_.clear();
    }

    ~Dispatcher() {
        wait();
    }
};

template <typename Task, typename... Args>
void dispatch(int affinity_id, int num_threads, Task task, Args... args) {
    Dispatcher dispatcher(affinity::cpu_assignment(affinity_id, num_threads), num_threads, task, args...);
    dispatcher.wait();
}

}  // namespace thread_coordination
