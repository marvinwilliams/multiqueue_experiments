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

template <bool Min, typename Key = unsigned long, typename T = Key>
class CAPQ {
    static_assert(std::is_same_v<Key, unsigned long> && std::is_same_v<T, unsigned long>,
                  "CA-PQ only supports unsigned long as key and value type");

   public:
    using key_type = unsigned long;
    using mapped_type = unsigned long;
    using value_type = std::pair<key_type, mapped_type>;

   private:
    struct Deleter {
        void operator()(::fpasl_catree_set* /*unused*/) {
            ::_destroy_gc_subsystem();
        }
    };
    alignas(64) std::unique_ptr<::fpasl_catree_set, Deleter> pq_;
    static constexpr key_type sentinel = std::numeric_limits<key_type>::max();

   public:
    using handle_type = SelfHandle<CAPQ>;
    using settings_type = util::EmptySettings;

    explicit CAPQ(int /*unused*/, std::size_t /*unused*/, settings_type const& /*unused*/) {
        ::_init_gc_subsystem();
        pq_.reset(::capq_new());
    }

    void push(value_type const& value) {
        ::capq_put_param(pq_.get(), Min ? value.first : sentinel - value.first - 1, value.second, true);
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

    static void write_human_readable(std::ostream& out) {
        out << "CA-PQ\n";
    }

    handle_type get_handle() {
        return handle_type{*this};
    }
};

}  // namespace wrapper::capq
