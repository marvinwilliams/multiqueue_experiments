#pragma once

// Adapted from klsm

#include "k_lsm/k_lsm.h"

#include "util.hpp"

#include <cxxopts.hpp>

#include <algorithm>
#include <iostream>
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
   private:
    using pq_type = ::kpq::k_lsm<key_type, mapped_type, relaxation>;
    static constexpr key_type sentinel_ = std::numeric_limits<key_type>::max();

   public:
    using handle_type = util::SelfHandle<KLsm>;
    using settings_type = util::EmptySettings;

   private:
    pq_type pq_{};

   public:
    explicit KLsm(int /*unused*/, std::size_t /*unused*/, settings_type const& /*unused*/) {};

    void push(value_type const& value) {
        pq_.insert(Min ? value.first : sentinel_ - value.first - 1, value.second);
    }

    std::optional<value_type> try_pop() {
        key_type key;
        mapped_type value;
        if (!pq_.delete_min(key, value)) {
            return std::nullopt;
        }
        return value_type{Min ? key : sentinel_ - key - 1, value};
    }

    static void write_human_readable(std::ostream& out) {
        out << "k-Lsm\n";
        out << "  k: " << relaxation << '\n';
    }

    handle_type get_handle() {
        return handle_type{*this};
    }
};

}  // namespace wrapper::klsm
