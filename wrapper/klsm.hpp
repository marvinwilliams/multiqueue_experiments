#pragma once

// Adapted from klsm

#include "wrapper/priority.hpp"

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

// The k-LSM can find remaining elements only from the thread that has the elements in its local buffer

template <typename Key, typename T, Priority P>
class KLSM {
   public:
    using key_type = Key;
    using mapped_type = T;
    using value_type = std::pair<key_type, mapped_type>;
#ifndef KLSM_K
    static constexpr std::size_t relaxation = 256;
#else
    static constexpr std::size_t relaxation = KLSM_K;
#endif
    struct config_type {};

   private:
    using pq_type = ::kpq::k_lsm<key_type, mapped_type, relaxation>;
    static constexpr key_type sentinel_ = std::numeric_limits<key_type>::max();

   public:
    class Handle {
        friend KLSM;
        pq_type* pq_;

        Handle(pq_type* pq) : pq_{pq} {
        }

       public:
        void push(value_type const& value) {
            pq_->insert(P == Priority::Min ? value.first : sentinel_ - value.first - 1, value.second);
        }

        std::optional<value_type> try_pop() {
            value_type retval;
            if (!pq_->delete_min(retval.first, retval.second)) {
                return std::nullopt;
            }
            if (P == Priority::Max) {
                retval.first = sentinel_ - retval.first - 1;
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

    KLSM(int /*num_threads*/, std::size_t /*initial_capacity*/, config_type const& /*options*/) : pq_(new pq_type()) {
    }

    Handle get_handle() {
        return Handle{pq_.get()};
    }

    std::ostream& describe(std::ostream& out) {
        out << "k-LSM (k: " << relaxation << ')';
        return out;
    }
};

}  // namespace wrapper
