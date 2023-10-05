#include "linden.hpp"

#include "wrapper/priority.hpp"

extern "C" {
#include "spraylist_linden/gc/gc.h"
#include "spraylist_linden/linden.h"
}
#undef min
#undef max

#include <cstddef>
#include <iostream>
#include <utility>

namespace wrapper {

struct detail::LindenBase::PQWrapper : ::pq_t {};

void detail::LindenBase::Deleter::operator()(PQWrapper* p) {
    // Avoid segfault
    ::insert(p, 1, 1);
    ::pq_destroy(p);
    ::_destroy_gc_subsystem();
}

detail::LindenBase::LindenBase(int /* num_threads */, std::size_t /*unused*/, config_type const& /*options*/) {
    _init_gc_subsystem();
    pq_.reset(static_cast<PQWrapper*>(pq_init(32)));
}

void detail::LindenBase::push(value_type const& value) {
    ::insert(pq_.get(), value.first, value.second);
}

auto detail::LindenBase::try_pop() -> std::optional<value_type> {
    value_type retval;
    retval.second = ::deletemin_key(pq_.get(), &retval.first);
    if (retval.first == sentinel) {
        return std::nullopt;
    }
    return retval;
}

}  // namespace wrapper
