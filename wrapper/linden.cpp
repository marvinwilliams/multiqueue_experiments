#include "linden.hpp"

#include <cstddef>
#include <utility>

extern "C" {
#include "spraylist_linden/gc/gc.h"
#include "spraylist_linden/linden.h"
}

namespace wrapper {

Linden::Linden(unsigned int /* num_threads */)
    : pq_((_init_gc_subsystem(), pq_init(32)), [](pq_t* pq) {
        // Avoid segfault
        ::insert(pq, 1, 1);
        pq_destroy(pq);
        _destroy_gc_subsystem();
      }) {}

Linden::Handle Linden::get_handle() {
  auto h = Handle{};
  h.pq_ = pq_.get();
  return h;
}

void Linden::Handle::push(value_type const& value) {
  ::insert(pq_, value.first + 1, value.second);
}

bool Linden::Handle::try_extract_top(value_type& retval) {
  unsigned long k_ret;
  retval.second = deletemin_key(pq_, &retval.first);
  --retval.first;
  return retval.first != sentinel_;
}

}  // namespace wrapper
