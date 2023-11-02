#pragma once

// Adapted from klsm

extern "C" {
#include "spraylist_linden/include/random.h"
#include "spraylist_linden/intset.h"
#include "ssalloc.h"
}

#undef min
#undef max

#include "cxxopts.hpp"

#include <cstdint>
#include <limits>
#include <memory>
#include <mutex>
#include <optional>
#include <ostream>
#include <utility>

namespace wrapper::spraylist {

template <bool Min>
class Spraylist {
   public:
    using key_type = unsigned long;
    using mapped_type = unsigned long;
    using value_type = std::pair<key_type, mapped_type>;

   private:
    struct PQWrapper;
    struct Deleter {
        void operator()(PQWrapper* p) {
            ::sl_set_delete(p);
        }
    };
    alignas(64) std::unique_ptr<PQWrapper, Deleter> pq_;
    int num_threads_;

    static constexpr key_type sentinel = std::numeric_limits<key_type>::max();

    struct ThreadDataWrapper : ::thread_data_t {};
    struct ThreadDataDeleter {
        void operator()(ThreadDataWrapper* p) {
            delete static_cast<thread_data_t*>(p);
        }
    };

    class Handle {
        friend Spraylist;
        Spraylist* pq_;
        std::unique_ptr<ThreadDataWrapper, ThreadDataDeleter> data_;

        explicit Handle(Spraylist* pq) : pq_{pq}, data_{pq_->new_thread_data()} {
        }

       public:
        void push(value_type const& value) {
            if constexpr (Min) {
                ::sl_add_val(pq_, value.first, value.second, TRANSACTIONAL);
            } else {
                ::sl_add_val(pq_, sentinel - value.first - 1, value.second, TRANSACTIONAL);
            }
        }

        std::optional<value_type> try_pop() {
            unsigned long key;
            unsigned long value;
            int ret;
            do {
                ret = ::spray_delete_min_key(pq_, &key, &value, data_.get());
            } while (ret == 0 && key != sentinel);
            if (key == sentinel) {
                return std::nullopt;
            }
            if constexpr (!Min) {
                key = sentinel - key - 1;
            }
            return value_type{key, value};
        };
    };

    [[nodiscard]] std::unique_ptr<ThreadDataWrapper, ThreadDataDeleter> new_thread_data() const {
        static std::mutex m;
        auto l = std::scoped_lock(m);
        ::ssalloc_init(num_threads_);
        seeds = ::seed_rand();
        auto thread_data =
            std::unique_ptr<ThreadDataWrapper, ThreadDataDeleter>(static_cast<ThreadDataWrapper*>(new thread_data_t));
        thread_data->seed = ::rand();
        thread_data->seed2 = ::rand();
        thread_data->nb_threads = num_threads_;
        return thread_data;
    }

   public:
    using handle_type = Handle;

    Spraylist(int num_threads, std::size_t initial_capacity) : num_threads_(num_threads) {
        ::ssalloc_init(num_threads_);
        *levelmax = floor_log_2(initial_capacity);
        pq_.reset(static_cast<PQWrapper*>(sl_set_new()));
    }

    Handle get_handle() {
        return Handle{this};
    }
};

template <bool Min>
using PQWrapper = Spraylist<Min>;

inline void add_options(cxxopts::Options& /*options*/) {
}

template <bool Min>
Spraylist<Min> create(int num_threads, std::size_t initial_capacity, cxxopts::ParseResult const& /*result*/) {
    return Spraylist<Min>{num_threads, initial_capacity};
}

template <bool Min>
std::ostream& describe(Spraylist<Min> const& /*unused*/, std::ostream& out) {
    out << "Spraylist";
    return out;
}

}  // namespace wrapper::spraylist
