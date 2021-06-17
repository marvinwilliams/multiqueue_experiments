#ifndef WRAPPER_LINDEN_HPP_INCLUDED
#define WRAPPER_LINDEN_HPP_INCLUDED

// Adapted from klsm

extern "C" {
#include "spraylist_linden/gc/gc.h"
#include "spraylist_linden/linden.h"
}

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>

namespace wrapper {

struct linden_pq_wrapper {
    pq_t* pq;
};

template <typename Key, typename T>
class linden {
    linden_pq_wrapper* pq_;

   public:
    static_assert(std::is_unsigned_v<Key> && std::is_unsigned_v<T>, "Only unsigned integers allowed");
    static_assert((std::numeric_limits<Key>::digits <= std::numeric_limits<unsigned long>::digits) &&
                      (std::numeric_limits<T>::digits <= std::numeric_limits<unsigned long>::digits),
                  "Type must not be larger than unsigned long");

    struct Handle {};
    static constexpr int DEFAULT_OFFSET = 32;

    linden(unsigned int num_threads = 0, int const max_offset = DEFAULT_OFFSET) {
        _init_gc_subsystem();
        pq_ = new linden_pq_wrapper;
        pq_->pq = pq_init(max_offset);
    }

    constexpr Handle get_handle(unsigned int = 0) {
        return Handle{};
    }

    ~linden() {
        // Avoid segfault
        push(Handle{}, {1u, 1u});
        pq_destroy(pq_->pq);
        delete pq_;
        _destroy_gc_subsystem();
    }

    void push(Handle, std::pair<Key, T> const& value) {
        ::insert(pq_->pq, static_cast<unsigned long>(value.first + 1), static_cast<unsigned long>(value.second));
    }

    bool extract_top(Handle, std::pair<Key, T>& retval) {
        unsigned long k_ret;
        retval.second = static_cast<T>(::deletemin_key(pq_->pq, &k_ret));
        retval.first = static_cast<Key>(k_ret - 1);
        return k_ret != -1;
    }

    static std::string description() {
        return "linden";
    }
};

}  // namespace wrapper

#endif
