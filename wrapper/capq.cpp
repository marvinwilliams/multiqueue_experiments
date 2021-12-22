#include "wrapper/capq.hpp"

extern "C" {
#include "capq/capq.h"
#include "capq/gc/gc.h"
}

__thread ptst_t* ptst;

namespace wrapper {

template <bool A, bool B, bool C>
struct Capq<A, B, C>::wrapper_type {
  char pad1[64 - sizeof(CAPQ*)];
  CAPQ* pq;
  char pad2[64];

  ~wrapper_type() {
    // freeing the capq segfaults
    // capq_delete(pq);
  }
};

template <bool A, bool B, bool C>
Capq<A, B, C>::Capq() : pq_(new wrapper_type) {
  _init_gc_subsystem();
  pq_->pq = capq_new();
}

template <bool A, bool B, bool C>
Capq<A, B, C>::~Capq() {
  pq_.reset();
  _destroy_gc_subsystem();
}

template <bool A, bool B, bool catree_adapt>
void Capq<A, B, catree_adapt>::push(value_type const& value) {
  capq_put_param(pq_->pq, value.first, value.second, catree_adapt);
}

template <bool remove_min_relax, bool put_relax, bool catree_adapt>
bool Capq<remove_min_relax, put_relax, catree_adapt>::try_extract_top(
    value_type& retval) {
  retval.second = capq_remove_min_param(pq_->pq, &retval.first, remove_min_relax,
                                      put_relax, catree_adapt);
  return retval.first != empty_key;
}

template class Capq<true, true, true>;
template class Capq<true, false, true>;
template class Capq<false, true, true>;
template class Capq<false, false, true>;

}  // namespace wrapper
