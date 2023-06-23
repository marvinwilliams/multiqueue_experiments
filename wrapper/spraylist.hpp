#pragma once

// Adapted from klsm

#include "cxxopts.hpp"

#include <optional>
#include <cstdint>
#include <limits>
#include <memory>
#include <ostream>
#include <utility>

using sl_intset_t = struct sl_intset;
using thread_data_t = struct thread_data;

namespace wrapper {

template <typename Key = unsigned long, typename T = unsigned long, bool Min = true>
class Spraylist {
    static_assert(std::is_same_v<Key, unsigned long> && std::is_same_v<T, unsigned long>,
                  "Only unsigned long as Key and T are supported");
};

template <bool Min>
class Spraylist<unsigned long, unsigned long, Min> {
    struct sl_intset_deleter {
        void operator()(sl_intset_t*);
    };
    struct thread_data_deleter {
        void operator()(thread_data_t*);
    };

   public:
    using key_type = unsigned long;
    using mapped_type = unsigned long;
    using value_type = std::pair<key_type, mapped_type>;

    struct Handle {
        std::unique_ptr<thread_data_t, thread_data_deleter> data_;
        sl_intset_t* pq_;

        void push(value_type const& value);
        std::optional<value_type> try_pop();
    };

    using handle_type = Handle;

    // Uses INT_MAX_32
    static constexpr key_type max_valid_key = std::numeric_limits<std::uint32_t>::max() - 1UL;
    static_assert(std::numeric_limits<unsigned long>::max() == -1UL);

   private:
    static constexpr key_type sentinel_ = std::numeric_limits<key_type>::max();

    alignas(64) std::unique_ptr<sl_intset_t, sl_intset_deleter> pq_;
    int num_threads_;

   public:
    static void add_options(cxxopts::Options& /*options*/) {
    }

    explicit Spraylist(int num_threads, std::size_t initial_capacity, cxxopts::ParseResult const& options);

    Handle get_handle();

    std::ostream& describe(std::ostream& out) {
        out << "Spraylist";
        return out;
    }
};

}  // namespace wrapper
