#pragma once

// Adapted from klsm

extern "C" {
#include "spraylist_linden/gc/gc.h"
#include "spraylist_linden/linden.h"
}
#undef min
#undef max

#include "cxxopts.hpp"

#include <limits>
#include <memory>
#include <optional>
#include <ostream>
#include <utility>

namespace wrapper::linden {

// The queue itself only supports
// keys >= 1, so one is added on each insert
template <bool Min>
class Linden {
   public:
    using key_type = unsigned long;
    using mapped_type = unsigned long;
    using value_type = std::pair<key_type, mapped_type>;
    struct config_type {};
    using handle_type = Linden&;

   private:
    struct Deleter {
        void operator()(::pq_t* p) {
            // Avoid segfault
            ::insert(p, 1, 1);
            ::pq_destroy(p);
            ::_destroy_gc_subsystem();
        };
    };

    alignas(64) std::unique_ptr<::pq_t, Deleter> pq_;

    static constexpr key_type sentinel = std::numeric_limits<key_type>::max();

   public:
    explicit Linden() {
        _init_gc_subsystem();
        pq_.reset(::pq_init(32));
    }

    void push(value_type const& value) {
        if constexpr (Min) {
            ::insert(pq_.get(), value.first + 1, value.second);
        } else {
            ::insert(pq_.get(), sentinel - value.first - 1, value.second);
        }
    }
    std::optional<value_type> try_pop() {
        unsigned long key;
        unsigned long value = ::deletemin_key(pq_.get(), &key);
        if (key == sentinel) {
            return std::nullopt;
        }
        if constexpr (Min) {
            --key;
        } else {
            key = sentinel - key - 1;
        }
        return value_type{key, value};
    }

    handle_type get_handle() {
        return *this;
    }
};

template <bool Min = true>
using PQWrapper = Linden<Min>;

inline void add_options(cxxopts::Options& /*options*/) {
}

template <bool Min = true>
Linden<Min> create(int /*num_threads*/, std::size_t /*initial_capacity*/, cxxopts::ParseResult const& /*result*/) {
    return Linden<Min>{};
}

template <typename PQ>
std::ostream& describe(PQ const& /*unused*/, std::ostream& out) {
    out << "Linden";
    return out;
}

}  // namespace wrapper::linden
