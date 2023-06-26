#include "wrapper/capq.hpp"
#include "capq.hpp"

extern "C" {
#include "capq/capq.h"
#include "capq/gc/gc.h"
}

#include <cstddef>

namespace wrapper {

void detail::CAPQ_deleter::operator()(capq_type* p) {
    _destroy_gc_subsystem();
}

template <bool Min>
CAPQ<unsigned long, unsigned long, Min>::CAPQ(int /*num_threads*/, std::size_t /*initial_capacity*/,
                                              config_type const& /*options*/) {
    _init_gc_subsystem();
    pq_.reset(capq_new());
}

template <>
void CAPQ<unsigned long, unsigned long, true>::Handle::push(value_type const& value) {
    capq_put_param(pq_, value.first, value.second, true);
}

template <>
void CAPQ<unsigned long, unsigned long, false>::Handle::push(value_type const& value) {
    capq_put_param(pq_, sentinel_ - 1 - value.first, value.second, true);
}

template <>
auto CAPQ<unsigned long, unsigned long, true>::Handle::try_pop() -> std::optional<value_type> {
    value_type retval;
    retval.second = capq_remove_min_param(pq_, &retval.first, true, true, true);
    if (retval.first == sentinel_) {
        return std::nullopt;
    }
    return retval;
}

template <>
auto CAPQ<unsigned long, unsigned long, false>::Handle::try_pop() -> std::optional<value_type> {
    value_type retval;
    retval.second = capq_remove_min_param(pq_, &retval.first, true, true, true);

    if (retval.first == sentinel_) {
        return std::nullopt;
    }
    retval.first = sentinel_ - 1 - retval.first;
    return retval;
}

template class CAPQ<unsigned long, unsigned long, true>;
template class CAPQ<unsigned long, unsigned long, false>;

}  // namespace wrapper
