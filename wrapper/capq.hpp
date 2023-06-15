#pragma once

// Adapted from klsm

#include "cxxopts.hpp"

#include <ios>
#include <limits>
#include <memory>
#include <ostream>
#include <utility>

using CAPQ = struct fpasl_catree_set;
struct CAPQ_deleter {
    void operator()(CAPQ*);
};

namespace wrapper {

// GC has MAX_THREADS = 128
// (unsigned long) -1 is signaling empty

template <typename Key = unsigned long, typename T = unsigned long, bool Min = true>
class CAPQ {
    static_assert(std::is_same_v<Key, unsigned long> && std::is_same_v<T, unsigned long>,
                  "Only unsigned long as Key and T are supported");
};

template <bool Min>
class CAPQ<unsigned long, unsigned long, Min> {
   public:
    using key_type = unsigned long;
    using mapped_type = unsigned long;
    using value_type = std::pair<key_type, mapped_type>;

    struct Handle {
        CAPQ* pq_;

        void push(value_type const& value);
        bool try_pop(value_type& retval);
    };

    using handle_type = Handle;

    static constexpr key_type min_valid_key = std::numeric_limits<key_type>::min();
    static constexpr key_type max_valid_key = std::numeric_limits<key_type>::max() - 1;

   private:
    static constexpr key_type sentinel_ = std::numeric_limits<key_type>::max();

    alignas(64) std::unique_ptr<CAPQ, CAPQ_deleter> pq_;

   public:
    static void add_options(cxxopts::Options& /*options*/) {
    }

    CAPQ(int num_threads, std::size_t initial_capacity, cxxopts::ParseResult const& options);

    Handle get_handle();

    std::ostream& describe(std::ostream& out) {
        out << "CA-PQ";
        return out;
    }
};

}  // namespace wrapper
