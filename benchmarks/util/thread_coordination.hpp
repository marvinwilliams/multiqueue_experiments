#pragma once

#include "barrier.hpp"
#include "threading.hpp"

#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <mutex>
#include <numeric>
#include <ostream>
#include <thread>
#include <utility>
#include <vector>

namespace thread_coordination {

namespace affinity {

namespace detail {
std::vector<std::array<std::size_t, 4>> get_cache_hierarchy() {
    std::vector<std::array<std::size_t, 4>> hierarchy(std::thread::hardware_concurrency());
    std::vector<std::size_t> l1_lookup;
    for (std::size_t i = 0; i < hierarchy.size(); ++i) {
        hierarchy[i][0] = i;
        auto ss = std::ifstream("/sys/devices/system/cpu/cpu" + std::to_string(i) + "/cache/index0/id");
        ss >> hierarchy[i][1];
        ss = std::ifstream("/sys/devices/system/cpu/cpu" + std::to_string(i) + "/cache/index3/id");
        ss >> hierarchy[i][2];
        ss = std::ifstream("/sys/devices/system/cpu/cpu" + std::to_string(i) + "/topology/physical_package_id");
        ss >> hierarchy[i][3];
    }
    for (std::size_t l = 0; l < 3; ++l) {
        std::vector<std::vector<std::size_t>> lookup;
        for (auto& h : hierarchy) {
            if (h[l + 1] >= lookup.size()) {
                lookup.resize(h[l + 1] + 1);
            }
            auto it = std::find(lookup[h[l + 1]].begin(), lookup[h[l + 1]].end(), h[l]);
            if (it == lookup[h[l + 1]].end()) {
                lookup[h[l + 1]].push_back(h[l]);
                it = lookup[h[l + 1]].end() - 1;
            }
            h[l] = std::size_t(std::distance(lookup[h[l + 1]].begin(), it));
        }
    }
    return hierarchy;
}
}  // namespace detail

struct None {
    threading::thread_config operator()(int /*unused*/) const {
        return {};
    }
};

struct ThreadId {
    std::size_t stride = 1;
    std::size_t offset = 0;

    threading::thread_config operator()(int id) const {
        threading::thread_config cfg;
        cfg.cpu_set.set(offset + static_cast<std::size_t>(id) * stride);
        return cfg;
    }
};

struct Same {
    std::size_t cpu = 1;

    threading::thread_config operator()(int /*unused*/) const {
        threading::thread_config cfg;
        cfg.cpu_set.set(cpu);
        return cfg;
    }
};

struct CloseCaches {
    std::vector<std::size_t> order;

    CloseCaches() {
        auto hierarchy = detail::get_cache_hierarchy();
        order.resize(hierarchy.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&hierarchy](std::size_t a, std::size_t b) {
            return std::tie(hierarchy[a][3], hierarchy[a][2], hierarchy[a][1], hierarchy[a][0]) <
                std::tie(hierarchy[b][3], hierarchy[b][2], hierarchy[b][1], hierarchy[b][0]);
        });
    }

    threading::thread_config operator()(int id) const {
        threading::thread_config cfg;
        cfg.cpu_set.set(order[std::size_t(id)]);
        return cfg;
    }
};

struct FarCaches {
    std::vector<std::size_t> order;

    FarCaches() {
        auto hierarchy = detail::get_cache_hierarchy();
        order.resize(hierarchy.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&hierarchy](std::size_t a, std::size_t b) {
            return std::tie(hierarchy[a][3], hierarchy[a][0], hierarchy[a][1], hierarchy[a][2]) <
                std::tie(hierarchy[b][3], hierarchy[b][0], hierarchy[b][1], hierarchy[b][2]);
        });
    }

    threading::thread_config operator()(int id) const {
        threading::thread_config cfg;
        cfg.cpu_set.set(order[std::size_t(id)]);
        return cfg;
    }
};

struct CloseL3FarL1 {
    std::vector<std::size_t> order;

    CloseL3FarL1() {
        auto hierarchy = detail::get_cache_hierarchy();
        order.resize(hierarchy.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&hierarchy](std::size_t a, std::size_t b) {
            return std::tie(hierarchy[a][3], hierarchy[a][2], hierarchy[a][0], hierarchy[a][1]) <
                std::tie(hierarchy[b][3], hierarchy[b][2], hierarchy[b][0], hierarchy[b][1]);
        });
    }

    threading::thread_config operator()(int id) const {
        threading::thread_config cfg;
        cfg.cpu_set.set(order[std::size_t(id)]);
        return cfg;
    }
};

struct FarL1CloseL3 {
    std::vector<std::size_t> order;

    FarL1CloseL3() : order(std::thread::hardware_concurrency()) {
        auto hierarchy = detail::get_cache_hierarchy();
        order.resize(hierarchy.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&hierarchy](std::size_t a, std::size_t b) {
            return std::tie(hierarchy[a][3], hierarchy[a][0], hierarchy[a][2], hierarchy[a][1]) <
                std::tie(hierarchy[b][3], hierarchy[b][0], hierarchy[b][2], hierarchy[b][1]);
        });
    }

    threading::thread_config operator()(int id) const {
        threading::thread_config cfg;
        cfg.cpu_set.set(order[std::size_t(id)]);
        return cfg;
    }
};

}  // namespace affinity

class Context {
    friend class Dispatcher;

    int id_;
    Barrier* barrier_;
    std::mutex* write_mutex_{};

    Context(int id, Barrier& barrier, std::mutex& m) : id_{id}, barrier_{&barrier}, write_mutex_{&m} {
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

    template <typename CompletionFunc>
    void synchronize(CompletionFunc&& f) const {
        barrier_->wait(std::forward<CompletionFunc>(f));
    }
};

class Dispatcher {
    std::vector<threading::pthread> threads_;
    Barrier barrier_;
    std::mutex write_mutex_;

   public:
    template <typename Affinity, typename Task, typename... Args>
    explicit Dispatcher(Affinity const& affinity, int num_threads, Task task, Args... args)
        : threads_(static_cast<std::size_t>(num_threads)), barrier_(num_threads) {
        for (int i = 0; i < num_threads; ++i) {
            threads_[static_cast<std::size_t>(i)] =
                threading::pthread(affinity(i), task, Context{i, barrier_, write_mutex_}, args...);
        }
    }

    template <typename Task, typename... Args>
    explicit Dispatcher(int num_threads, Task task, Args&&... args)
        : Dispatcher(affinity::FarL1CloseL3{}, num_threads, task, std::forward<Args>(args)...) {
    }

    void wait() {
        for (auto& t : threads_) {
            t.join();
        }
    }
};

}  // namespace thread_coordination
