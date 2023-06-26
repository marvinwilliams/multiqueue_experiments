#pragma once

// Adapted from klsm

#include "cxxopts.hpp"

#include <limits>
#include <memory>
#include <optional>
#include <ostream>
#include <utility>

namespace wrapper {

template <typename Key = unsigned long, typename T = unsigned long, bool Min = true>
class Linden {
    static_assert(std::is_same_v<Key, unsigned long> && std::is_same_v<T, unsigned long>,
                  "Only unsigned long as Key and T are supported");
};

template <bool Min>
class Linden<unsigned long, unsigned long, Min> {
    struct linden_pq_t;
    struct pq_deleter {
        void operator()(linden_pq_t*);
    };

   public:
    using key_type = unsigned long;
    using mapped_type = unsigned long;
    using value_type = std::pair<key_type, mapped_type>;
    struct config_type {};

    struct Handle {
        linden_pq_t* pq_;

        void push(value_type const& value);
        std::optional<value_type> try_pop();
    };

    using handle_type = Handle;

    // The queue itself only supports
    // keys >= 1, so one is added on each insert
    static constexpr key_type max_valid_key = std::numeric_limits<key_type>::max() >> 1;

   private:
    static constexpr key_type sentinel_ = std::numeric_limits<key_type>::max();

    struct wrapper_type;

    alignas(64) std::unique_ptr<linden_pq_t, pq_deleter> pq_;

   public:
    static void add_options(cxxopts::Options& /*options*/, config_type& /*config*/) {
    }

    Linden(int num_threads, std::size_t initial_capacity, config_type const& options);

    Handle get_handle();

    std::ostream& describe(std::ostream& out) {
        out << "Linden";
        return out;
    }
};

}  // namespace wrapper
