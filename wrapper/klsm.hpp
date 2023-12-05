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

namespace wrapper::klsm {

// The k-LSM can find remaining elements only from the thread that has the elements in its local buffer

template <bool Min, typename Key = unsigned long, typename T = unsigned long>
class KLsm {
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
    using handle_type = KLsm&;

   private:
    using pq_type = ::kpq::k_lsm<key_type, mapped_type, relaxation>;
    static constexpr key_type sentinel_ = std::numeric_limits<key_type>::max();

   public:
   private:
    pq_type pq_{};

   public:
    explicit KLsm() = default;

    void push(value_type const& value) {
        pq_.insert(Min ? value.first : sentinel_ - value.first - 1, value.second);
    }

    std::optional<value_type> try_pop() {
        key_type key;
        mapped_type value;
        if (!pq_.delete_min(key, value)) {
            return std::nullopt;
        }
        if (!Min) {
            key = sentinel_ - key - 1;
        }
        return value_type{key, value};
    }

    handle_type get_handle() {
        return *this;
    }
};

template <bool Min = true, typename Key = unsigned long, typename Value = std::pair<unsigned long, unsigned long>>
using PQWrapper = KLsm<Min, Key, T>;

inline void add_options(cxxopts::Options& /*options*/) {
}

template <bool Min = true, typename Key = unsigned long, typename T = unsigned long>
KLsm<Min, Key, T> create(int /*num_threads*/, std::size_t /*initial_capacity*/,
                         cxxopts::ParseResult const& /*result*/) {
    return KLsm<Min, Key, T>{};
}

template <typename PQ>
std::ostream& describe(PQ const& /*unused*/, std::ostream& out) {
    out << "k-LSM (k: " << PQ::relaxation << ')';
    return out;
}

}  // namespace wrapper::klsm
