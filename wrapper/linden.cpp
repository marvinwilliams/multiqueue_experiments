#include "linden.hpp"

#include <cstddef>
#include <utility>

extern "C" {
#include "spraylist_linden/gc/gc.h"
#include "spraylist_linden/linden.h"
}

namespace wrapper {

struct Linden::wrapper_type {
  pq_t* pq;
  ~wrapper_type() { pq_destroy(pq); }
};

Linden::Linden() : pq_(new wrapper_type) {
  _init_gc_subsystem();
  pq_->pq = pq_init(32);
}

Linden::~Linden() {
  // Avoid segfault
  ::insert(pq_->pq, 1, 1);
  pq_.reset();
  _destroy_gc_subsystem();
}

void Linden::push(Handle&, value_type value) {
  ::insert(pq_->pq, value.key + 1, value.data);
}

bool Linden::try_delete_min(Handle&, value_type& retval) {
  unsigned long k_ret;
  retval.data = deletemin_key(pq_->pq, &retval.key);
  --retval.key;
  return retval.key != empty_key;
}

}  // namespace wrapper
