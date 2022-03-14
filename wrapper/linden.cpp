#include "linden.hpp"

extern "C" {
#include "spraylist_linden/gc/gc.h"
#include "spraylist_linden/linden.h"
}
#undef min
#undef max

#include <cstddef>
#include <utility>

struct linden_pq_t : ::pq_t {};
void linden_pq_deleter::operator()(linden_pq_t* p) {
  // Avoid segfault
  ::insert(p, 1, 1);
  pq_destroy(p);
  _destroy_gc_subsystem();
}

namespace wrapper {

Linden::Linden(std::size_t /* capacity */, unsigned int /* num_threads */) {
  _init_gc_subsystem();
  pq_.reset(static_cast<linden_pq_t*>(pq_init(32)));
}

void Linden::Handle::push(value_type const& value) {
  ::insert(pq_, value.first + 1, value.second);
}

bool Linden::Handle::try_extract_top(value_type& retval) {
  unsigned long k_ret;
  retval.second = ::deletemin_key(pq_, &retval.first);
  --retval.first;
  return retval.first != sentinel_;
}

Linden::Handle Linden::get_handle() {
  auto h = Handle{};
  h.pq_ = pq_.get();
  return h;
}

void Linden::push(value_type const& value) {
  ::insert(pq_.get(), value.first + 1, value.second);
}

bool Linden::try_extract_top(value_type& retval) {
  unsigned long k_ret;
  retval.second = ::deletemin_key(pq_.get(), &retval.first);
  --retval.first;
  return retval.first != sentinel_;
}

}  // namespace wrapper
