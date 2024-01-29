#pragma once

// Adapted from klsm

extern "C" {
#include "spraylist_linden/gc/gc.h"
#include "spraylist_linden/linden.h"
}
#undef min
#undef max

#include <cxxopts.hpp>

#include <limits>
#include <memory>
#include <optional>
#include <ostream>
#include <utility>

namespace wrapper::linden {

// The queue itself only supports
// keys >= 1, so one is added on each insert
//
// Each key can only be inserted once
template <bool Min, typename Key = unsigned long, typename T = Key>
class Linden {
    static_assert(std::is_same_v<Key, unsigned long> && std::is_same_v<T, unsigned long>,
                  "Linden only supports unsigned long as key and value type");

   public:
    using key_type = unsigned long;
    using mapped_type = unsigned long;
    using value_type = std::pair<key_type, mapped_type>;

   private:
    static constexpr key_type sentinel = std::numeric_limits<key_type>::max();

    struct Deleter {
        void operator()(::pq_t* p) {
            // Avoid segfault
            ::insert(p, 1, 1);
            ::pq_destroy(p);
            ::_destroy_gc_subsystem();
        };
    };

    alignas(64) std::unique_ptr<::pq_t, Deleter> pq_;

   public:
    using handle_type = util::SelfHandle<Linden>;
    using settings_type = util::EmptySettings;
    explicit Linden(int /*unused*/, std::size_t /*unused*/, settings_type const& /*unused*/) {
        _init_gc_subsystem();
        pq_.reset(::pq_init(32));
    }

    void push(value_type const& value) {
        ::insert(pq_.get(), Min ? value.first + 1 : sentinel - value.first - 1, value.second);
    }

    std::optional<value_type> try_pop() {
        unsigned long key;
        unsigned long value = ::deletemin_key(pq_.get(), &key);
        if (key == sentinel) {
            return std::nullopt;
        }
        return value_type{Min ? key - 1 : sentinel - key - 1, value};
    }

    static void write_human_readable(std::ostream& out) {
        out << "Linden\n";
    }

    handle_type get_handle() {
        return handle_type{pq_.get()};
    }
};

}  // namespace wrapper::linden
