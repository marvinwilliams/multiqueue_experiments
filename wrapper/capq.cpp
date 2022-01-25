#include "wrapper/capq.hpp"

extern "C" {
#include "capq/capq.h"
#include "capq/gc/gc.h"
}

__thread ptst_t* ptst;

namespace wrapper {

template <bool remove_min_relax, bool put_relax, bool catree_adapt>
Capq<remove_min_relax, put_relax, catree_adapt>::Capq(
    unsigned int /* num_threads */)
    : pq_((_init_gc_subsystem(), capq_new()), [](CAPQ* pq) {
        capq_delete(pq);
        _destroy_gc_subsystem();
      }) {}

template <bool remove_min_relax, bool put_relax, bool catree_adapt>
typename Capq<remove_min_relax, put_relax, catree_adapt>::Handle
Capq<remove_min_relax, put_relax, catree_adapt>::get_handle() {
  auto h = Handle{};
  h.pq_ = pq_.get();
  return h;
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

template class Capq<true, true, true>;
template class Capq<true, false, true>;
template class Capq<false, true, true>;
template class Capq<false, false, true>;

}  // namespace wrapper
