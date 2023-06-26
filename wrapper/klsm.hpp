#pragma once

// Adapted from klsm

#include "k_lsm/k_lsm.h"

#include "cxxopts.hpp"

#include <algorithm>
#include <cstdio>
#include <limits>
#include <memory>
#include <optional>
#include <ostream>
#include <utility>

namespace wrapper {

// No known limitations

template <typename KeyType, typename T, bool Min, int Relaxation>
class KLsm {
   public:
    using key_type = KeyType;
    using mapped_type = T;
    using value_type = std::pair<key_type, mapped_type>;
    struct config_type {};

   private:
    using pq_type = kpq::k_lsm<key_type, mapped_type, Relaxation>;

   public:
    class Handle {
        pq_type& pq_;

       public:
        Handle(pq_type& pq) : pq_{pq} {
        }

        void push(value_type const& value) {
            pq_.insert(Min ? value.first : std::numeric_limits<key_type>::max() - value.first - 1, value.second);
        }

        std::optional<value_type> try_pop() {
            value_type retval;
            if (!pq_.delete_min(retval.first, retval.second)) {
                return std::nullopt;
            }
            if (!Min) {
                retval.first = std::numeric_limits<key_type>::max() - retval.first - 1;
            }
            return retval;
        }
    };

    using handle_type = Handle;

   private:
    alignas(64) std::unique_ptr<pq_type> pq_;

   public:
    static void add_options(cxxopts::Options& /*options*/, config_type& /*config*/) {
    }

    KLsm(int /*num_threads*/, std::size_t /*initial_capacity*/, config_type const& /*options*/)
        : pq_(new pq_type()) {
    }

    Handle get_handle() {
        return Handle{*pq_.get()};
    }

    std::ostream& describe(std::ostream& out) {
        out << "k-LSM (k: " << Relaxation << ')';
        return out;
    }
};

}  // namespace wrapper
