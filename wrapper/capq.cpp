#include "wrapper/capq.hpp"
#include "capq.hpp"

extern "C" {
#include "capq/capq.h"
#include "capq/gc/gc.h"
}

__thread ptst_t* ptst;

void CAPQ_deleter::operator()(CAPQ* p) {
  // Segfaults
  /* capq_delete(p); */
  _destroy_gc_subsystem();
}

namespace wrapper {

template <bool remove_min_relax, bool put_relax, bool catree_adapt>
Capq<remove_min_relax, put_relax, catree_adapt>::Capq(
    std::size_t /* capacity */, unsigned int /* num_threads */) {
  _init_gc_subsystem();
  pq_.reset(capq_new());
}

template <bool remove_min_relax, bool put_relax, bool catree_adapt>
typename Capq<remove_min_relax, put_relax, catree_adapt>::Handle
Capq<remove_min_relax, put_relax, catree_adapt>::get_handle() {
  return Handle{pq_.get()};
}

template <bool remove_min_relax, bool put_relax, bool catree_adapt>
void Capq<remove_min_relax, put_relax, catree_adapt>::Handle::push(
    value_type const& value) {
  capq_put_param(pq_, value.first, value.second, catree_adapt);
}

template <bool remove_min_relax, bool put_relax, bool catree_adapt>
bool Capq<remove_min_relax, put_relax, catree_adapt>::Handle::try_extract_top(
    value_type& retval) {
  retval.second = capq_remove_min_param(pq_, &retval.first, remove_min_relax,
                                        put_relax, catree_adapt);
  return retval.first != sentinel_;
}

template <bool remove_min_relax, bool put_relax, bool catree_adapt>
void Capq<remove_min_relax, put_relax, catree_adapt>::push(
    value_type const& value) {
  capq_put_param(pq_.get(), value.first, value.second, catree_adapt);
}

template <bool remove_min_relax, bool put_relax, bool catree_adapt>
bool Capq<remove_min_relax, put_relax, catree_adapt>::try_extract_top(
    value_type& retval) {
  retval.second = capq_remove_min_param(pq_.get(), &retval.first, remove_min_relax,
                                        put_relax, catree_adapt);
  return retval.first != sentinel_;
}

template class Capq<true, true, true>;
template class Capq<true, false, true>;
template class Capq<false, true, true>;
template class Capq<false, false, true>;

}  // namespace wrapper
