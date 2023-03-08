#pragma once

// Adapted from klsm

#include "k_lsm/k_lsm.h"

#include <algorithm>
#include <cstdio>
#include <limits>
#include <memory>
#include <ostream>
#include <utility>

namespace wrapper {

// No known limitations

template <typename KeyType, typename T, int Relaxation>
class Klsm {
   public:
    using key_type = KeyType;
    using mapped_type = T;
    using value_type = std::pair<key_type, mapped_type>;

   private:
    using pq_type = kpq::k_lsm<key_type, mapped_type, Relaxation>;

   public:
    class Handle {
        pq_type& pq_;

       public:
        Handle(pq_type& pq) : pq_{pq} {
        }

        void push(value_type const& value) {
            pq_.insert(value.first, value.second);
        }
        bool try_pop(value_type& retval) {
            return pq_.delete_min(retval.first, retval.second);
        }
    };

    using handle_type = Handle;

   private:
    alignas(64) std::unique_ptr<pq_type> pq_;

   public:
    Klsm(int num_threads, std::size_t initial_capacity);

    Handle get_handle(int id);

    static std::ostream& describe(std::ostream& out) {
        out << "klsm\n";
        out << "Relaxation: " << Relaxation << '\n';
        return out;
    }
};

template <typename KeyType, typename T, int Relaxation>
Klsm<KeyType, T, Relaxation>::Klsm(int /*unused*/, std::size_t /*unused*/) : pq_(new pq_type()) {
}

template <typename KeyType, typename T, int Relaxation>
typename Klsm<KeyType, T, Relaxation>::Handle Klsm<KeyType, T, Relaxation>::get_handle(int /*unused*/) {
    return Handle{*pq_.get()};
}
}  // namespace wrapper
