#pragma once

// Adapted from klsm

#include "k_lsm/k_lsm.h"

#include "cxxopts.hpp"

#include <algorithm>
#include <cstdio>
#include <limits>
#include <memory>
#include <ostream>
#include <utility>

namespace wrapper {

// No known limitations

template <typename KeyType, typename T, bool Min>
class KLsm {
   public:
    using key_type = KeyType;
    using mapped_type = T;
    using value_type = std::pair<key_type, mapped_type>;

   private:
#ifdef KLSM_K
    static constexpr auto relaxation = KLSM_K;
#else
    static constexpr auto relaxation = 256;
#endif
    using pq_type = kpq::k_lsm<key_type, mapped_type, relaxation>;

   public:
    class Handle {
        pq_type& pq_;

       public:
        Handle(pq_type& pq) : pq_{pq} {
        }

        void push(value_type const& value) {
            pq_.insert(Min ? value.first : std::numeric_limits<key_type>::max - value.first - 1, value.second);
        }

        bool try_pop(value_type& retval) {
            if (Min) {
                return pq_.delete_min(retval.first, retval.second);
            } else {
                if (!pq_.delete_min(retval.first, retval.second)) {
                    return false;
                }
                retval.first = std::numeric_limits<key_type>::max - retval.first - 1;
                return true;
            }
        }
    };

    using handle_type = Handle;

   private:
    alignas(64) std::unique_ptr<pq_type> pq_;

   public:
    static void add_options(cxxopts::Options& /*options*/) {
    }

    KLsm(int /*num_threads*/, std::size_t /*initial_capacity*/, cxxopts::ParseResult const& /*options*/)
        : pq_(new pq_type()) {
    }

    Handle get_handle() {
        return Handle{*pq_.get()};
    }

    std::ostream& describe(std::ostream& out) {
        out << "k-LSM (k: " << relaxation << ')';
        return out;
    }
};

}  // namespace wrapper
