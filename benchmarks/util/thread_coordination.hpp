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

struct automatic {
    threading::thread_config operator()(int /*unused*/) const {
        return {};
    }
};

struct distinct_cpu {
    std::size_t stride = 1;
    std::size_t offset = 0;

    threading::thread_config operator()(int id) const {
        threading::thread_config cfg;
        cfg.cpu_set.set(offset + static_cast<std::size_t>(id) * stride);
        return cfg;
    }
};

struct same_cpu {
    std::size_t cpu = 1;

    threading::thread_config operator()(int /*unused*/) const {
        threading::thread_config cfg;
        cfg.cpu_set.set(cpu);
        return cfg;
    }
};

struct same_caches {
    std::vector<std::size_t> order;

    same_caches() : order(std::thread::hardware_concurrency()) {
        std::vector<std::array<std::size_t, 4>> caches(order.size());
        for (std::size_t i = 0; i < caches.size(); ++i) {
            caches[i][0] = i;
            for (std::size_t l = 1; l < 4; ++l) {
                auto ss = std::ifstream("/sys/devices/system/cpu/cpu" + std::to_string(i) + "/cache/index" +
                                        std::to_string(l == 1 ? 0 : l) + "/id");
                ss >> caches[i][l];
            }
        }
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&caches](std::size_t a, std::size_t b) {
            return std::lexicographical_compare(caches[a].rbegin(), caches[a].rend(), caches[b].rbegin(), caches[b].rend());
        });
    }

    threading::thread_config operator()(int id) const {
        threading::thread_config cfg;
        cfg.cpu_set.set(order[std::size_t(id)]);
        return cfg;
    }
};

struct distinct_caches {
    std::vector<std::size_t> order;

    distinct_caches() : order(std::thread::hardware_concurrency()) {
        std::vector<std::array<std::size_t, 5>> caches(order.size());
        for (std::size_t i = 0; i < caches.size(); ++i) {
            caches[i][0] = i;
            for (std::size_t l = 1; l < 4; ++l) {
                auto ss = std::ifstream("/sys/devices/system/cpu/cpu" + std::to_string(i) + "/cache/index" +
                                        std::to_string(l == 1 ? 0 : l) + "/id");
                ss >> caches[i][l];
            }
        }
        for (std::size_t l = 0; l < 4; ++l) {
            std::size_t num_higher_caches =
                (*std::max_element(caches.begin(), caches.end(),
                                   [l](const auto& a, const auto& b) { return a[l + 1] < b[l + 1]; }))[l + 1] +
                1;
            std::vector<std::vector<std::size_t>> ids(num_higher_caches);
            for (auto& cache : caches) {
                auto it = std::find(ids[cache[l + 1]].begin(), ids[cache[l + 1]].end(), cache[l]);
                auto new_id = std::size_t(std::distance(ids[cache[l + 1]].begin(), it));
                if (it == ids[cache[l + 1]].end()) {
                    ids[cache[l + 1]].push_back(cache[l]);
                }
                cache[l] = new_id;
            }
        }
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&caches](std::size_t a, std::size_t b) {
            return std::lexicographical_compare(caches[a].begin(), caches[a].end(), caches[b].begin(), caches[b].end());
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
    explicit Dispatcher(Affinity affinity, int num_threads, Task task, Args... args)
        : threads_(static_cast<std::size_t>(num_threads)), barrier_(num_threads) {
        for (int i = 0; i < num_threads; ++i) {
            threads_[static_cast<std::size_t>(i)] =
                threading::pthread(affinity(i), task, Context{i, barrier_, write_mutex_}, args...);
        }
    }

    template <typename Task, typename... Args>
    explicit Dispatcher(int num_threads, Task task, Args&&... args)
        : Dispatcher(affinity::same_caches{}, num_threads, task, std::forward<Args>(args)...) {
    }

    void wait() {
        for (auto& t : threads_) {
            t.join();
        }
    }
};

}  // namespace thread_coordination
