#include "wrapper/capq.hpp"
#include "capq.hpp"

extern "C" {
#include "capq/capq.h"
#include "capq/gc/gc.h"
}

#include "cxxopts.hpp"

#include <cstddef>

void CAPQ_deleter::operator()(CAPQ* p) {
    _destroy_gc_subsystem();
}

namespace wrapper {

template <bool Min, bool remove_min_relax, bool put_relax, bool catree_adapt>
CAPQ<unsigned long, unsigned long, Min, remove_min_relax, put_relax, catree_adapt>::CAPQ(
    int /*num_threads*/, std::size_t /*initial_capacity*/, cxxopts::ParseResult const& /*options*/) {
    _init_gc_subsystem();
    pq_.reset(capq_new());
}

template <bool Min, bool remove_min_relax, bool put_relax, bool catree_adapt>
auto CAPQ<unsigned long, unsigned long, Min, remove_min_relax, put_relax, catree_adapt>::get_handle() {
    return Handle{pq_.get()};
}

template <bool Min, bool remove_min_relax, bool put_relax, bool catree_adapt>
void CAPQ<unsigned long, unsigned long, Min, remove_min_relax, put_relax, catree_adapt>::Handle::push(
    value_type const& value) {
    capq_put_param(pq_, value.first, value.second, catree_adapt);
}

template <bool remove_min_relax, bool put_relax, bool catree_adapt>
void CAPQ<unsigned long, unsigned long, false, remove_min_relax, put_relax, catree_adapt>::Handle::push(
    value_type const& value) {
    capq_put_param(pq_, sentinel_ - 1 - value.first, value.second, catree_adapt);
}

template <bool Min, bool remove_min_relax, bool put_relax, bool catree_adapt>
auto CAPQ<unsigned long, unsigned long, Min, remove_min_relax, put_relax, catree_adapt>::Handle::try_pop()
    -> std::optional<value_type> {
    value_type retval;
    retval.second = capq_remove_min_param(pq_, &retval.first, remove_min_relax, put_relax, catree_adapt);
    if (retval.first == sentinel_) {
        return std::nullopt;
    }
    return retval;
}

template <bool remove_min_relax, bool put_relax, bool catree_adapt>
auto CAPQ<unsigned long, unsigned long, false, remove_min_relax, put_relax, catree_adapt>::Handle::try_pop()
    -> std::optional<value_type> {
    value_type retval;
    retval.second = capq_remove_min_param(pq_, &retval.first, remove_min_relax, put_relax, catree_adapt);

    if (retval.first == sentinel_) {
        return std::nullopt;
    }
    retval.first = sentinel_ - 1 - retval.first;
    return retval;
}

template class CAPQ<unsigned long, unsigned long, true>;
template class CAPQ<unsigned long, unsigned long, false>;

}  // namespace wrapper
