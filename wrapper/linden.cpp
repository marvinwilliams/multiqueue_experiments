#include "linden.hpp"

extern "C" {
#include "spraylist_linden/gc/gc.h"
#include "spraylist_linden/linden.h"
}
#undef min
#undef max

#include "cxxopts.hpp"

#include <cstddef>
#include <iostream>
#include <utility>

namespace wrapper {

template <bool Min>
struct Linden<unsigned long, unsigned long, Min>::linden_pq_t : pq_t {};

template <bool Min>
void Linden<unsigned long, unsigned long, Min>::pq_deleter::operator()(linden_pq_t* p) {
    // Avoid segfault
    ::insert(static_cast<pq_t*>(p), 1, 1);
    pq_destroy(static_cast<pq_t*>(p));
    _destroy_gc_subsystem();
}

template <bool Min>
Linden<unsigned long, unsigned long, Min>::Linden(int /* num_threads */, std::size_t /*unused*/,
                                                  cxxopts::ParseResult const& /*options*/) {
    _init_gc_subsystem();
    pq_.reset(static_cast<linden_pq_t*>(pq_init(32)));
}

template <>
void Linden<unsigned long, unsigend long, true>::Handle::push(value_type const& value) const {
    ::insert(pq_, value.first + 1, value.second);
}

template <>
void Linden<unsigned long, unsigend long, false>::Handle::push(value_type const& value) const {
    ::insert(pq_, max_valid_key - value.first + 1, value.second);
}

template <>
bool Linden<unsigned long, unsigned long, true>::Handle::try_pop(value_type& retval) const {
    retval.second = ::deletemin_key(pq_, &retval.first);
    if (retval.first == sentinel_) {
        return false;
    }
    --retval.first;
    return true;
}

template <>
bool Linden<unsigned long, unsigned long, false>::Handle::try_pop(value_type& retval) const {
    retval.second = ::deletemin_key(pq_, &retval.first);
    if (retval.first == sentinel_) {
        return false;
    }
    retval.first = max_valid_key - value.first + 1;
    return true;
}

template <bool Min>
auto Linden<unsigned long, unsigned long, Min>::get_handle() {
    auto h = Handle{};
    h.pq_ = pq_.get();
    return h;
}

template class Linden<unsigned long, unsigned long, true>;
template class Linden<unsigned long, unsigned long, false>;

}  // namespace wrapper
