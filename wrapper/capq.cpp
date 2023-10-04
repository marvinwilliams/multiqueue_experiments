#include "wrapper/capq.hpp"

#include "wrapper/priority.hpp"

extern "C" {
#include "capq/capq.h"
#include "capq/gc/gc.h"
}

#include <cstddef>

namespace wrapper {

void detail::CAPQBase::Deleter::operator()(capq_type* p) {
    ::_destroy_gc_subsystem();
}

detail::CAPQBase::CAPQBase(int /*num_threads*/, std::size_t /*initial_capacity*/,
                                              config_type const& /*options*/) {
    ::_init_gc_subsystem();
    pq_.reset(capq_new());
}

void detail::CAPQBase::push(value_type const& value) {
    ::capq_put_param(pq_.get(), value.first, value.second, true);
}

auto detail::CAPQBase::try_pop() -> std::optional<value_type> {
    value_type retval;
    retval.second = ::capq_remove_min_param(pq_.get(), &retval.first, true, true, true);
    if (retval.first == sentinel) {
        return std::nullopt;
    }
    return retval;
}

}  // namespace wrapper
