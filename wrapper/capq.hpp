#pragma once

// Adapted from klsm

#include "cxxopts.hpp"

#include <ios>
#include <limits>
#include <memory>
#include <optional>
#include <ostream>
#include <utility>

extern "C" {
#include "capq/capq.h"
#include "capq/gc/gc.h"
}

namespace wrapper::capq {

// GC has MAX_THREADS = 128
// (unsigned long) -1 is signaling empty

template <bool Min = true>
class CAPQ {
   public:
    using key_type = unsigned long;
    using mapped_type = unsigned long;
    using value_type = std::pair<key_type, mapped_type>;
    using handle_type = CAPQ&;

   private:
    struct Deleter {
        void operator()(::fpasl_catree_set* /*unused*/) {
            ::_destroy_gc_subsystem();
        }
    };
    alignas(64) std::unique_ptr<::fpasl_catree_set, Deleter> pq_;

    static constexpr key_type sentinel = std::numeric_limits<key_type>::max();

   public:
    explicit CAPQ() {
        ::_init_gc_subsystem();
        pq_.reset(::capq_new());
    }

    void push(value_type const& value) {
        if constexpr (Min) {
            ::capq_put_param(pq_.get(), value.first, value.second, true);
        } else {
            ::capq_put_param(pq_.get(), sentinel - value.first - 1, value.second, true);
        }
    }
    std::optional<value_type> try_pop() {
        unsigned long key;
        unsigned long value = ::capq_remove_min_param(pq_.get(), &key, true, true, true);
        if (key == sentinel) {
            return std::nullopt;
        }
        if constexpr (!Min) {
            key = sentinel - key - 1;
        }
        return value_type{key, value};
    }

    handle_type get_handle() {
        return *this;
    }
};

template <bool Min = true>
using PQWrapper = CAPQ<Min>;

inline void add_options(cxxopts::Options& /*options*/) {
}

template <bool Min = true>
CAPQ<Min> create(int /*num_threads*/, std::size_t /*initial_capacity*/, cxxopts::ParseResult const& /*result*/) {
    return CAPQ<Min>{};
}

template <bool Min>
std::ostream& describe(CAPQ<Min> const& /*unused*/, std::ostream& out) {
    out << "CA-PQ";
    return out;
}

}  // namespace wrapper::capq
