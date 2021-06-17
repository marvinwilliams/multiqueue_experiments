#ifndef WRAPPER_CAPQ_HPP_INCLUDED
#define WRAPPER_CAPQ_HPP_INCLUDED

// Adapted from klsm

extern "C" {
#include "capq/capq.h"
#include "capq/gc/gc.h"
}

#include <cstddef>
#include <cstdint>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <type_traits>

namespace wrapper {

struct capq_wrapper_type {
    char pad1[64 - sizeof(CAPQ*)];
    CAPQ* pq;
    char pad2[64];
};

template <typename Key, typename T, bool remove_min_relax = true, bool put_relax = true, bool catree_adapt = true>
class capq {
    capq_wrapper_type* pq_;

   public:
    static_assert(std::is_unsigned<Key>::value && std::is_unsigned<T>::value, "Only unsigned integers allowed");
    static_assert((std::numeric_limits<Key>::digits <= std::numeric_limits<unsigned long>::digits) &&
                      (std::numeric_limits<T>::digits <= std::numeric_limits<unsigned long>::digits),
                  "Type must not be larger than unsigned long");

    struct Handle {};

    capq() {
        _init_gc_subsystem();
        /* init_thread(1); */
        pq_ = new capq_wrapper_type;
        pq_->pq = capq_new();
    }

    explicit capq(unsigned int) : capq() {
    }

    constexpr Handle get_handle(unsigned int) {
        return Handle{};
    }

    void push(Handle, std::pair<Key, T> const& value) {
        capq_put_param(pq_->pq, static_cast<unsigned long>(value.first), static_cast<unsigned long>(value.second),
                       catree_adapt);
    }

    bool extract_top(Handle, std::pair<Key, T>& retval) {
        unsigned long key_write_back;
        retval.second =
            static_cast<T>(capq_remove_min_param(pq_->pq, &key_write_back, remove_min_relax, put_relax, catree_adapt));
        retval.first = static_cast<Key>(key_write_back);
        return key_write_back != std::numeric_limits<unsigned long>::max();
    }

    static std::string description() {
        std::stringstream ss;
        ss << "capq\n";
        ss << "Remove min relax: " << std::boolalpha << remove_min_relax << '\n';
        ss << "Put relax" << std::boolalpha << put_relax << '\n';
        ss << "Catree adapt" << std::boolalpha << catree_adapt << '\n';
        return ss.str();
    }
};

}  // namespace wrapper

#endif
