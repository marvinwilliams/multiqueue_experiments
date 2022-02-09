// Adapted from klsm

#include "wrapper/spraylist.hpp"

extern "C" {
#include "spraylist_linden/include/random.h"
#include "spraylist_linden/intset.h"
#include "ssalloc.h"
}

#undef min
#undef max

#include <cstddef>
#include <iostream>
#include <memory>

__thread unsigned long* seeds;

void sl_intset_deleter::operator()(sl_intset_t* p) { sl_set_delete(p); }
void thread_data_deleter::operator()(thread_data_t* p) { delete p; }

namespace wrapper {

Spraylist::Spraylist(std::size_t /* capacity */, unsigned int num_threads) : num_threads_(num_threads) {
  ssalloc_init(num_threads_);
  *levelmax = floor_log_2(1'000'000);
  pq_.reset(sl_set_new());
}

Spraylist::Handle Spraylist::get_handle() {
  ssalloc_init(num_threads_);
  seeds = seed_rand();
  Handle handle{};
  handle.data_.reset(new thread_data_t());
  handle.data_->seed = rand();
  handle.data_->seed2 = rand();
  handle.data_->nb_threads = num_threads_;
  handle.pq_ = pq_.get();
  return handle;
}

void Spraylist::Handle::push(value_type const& value) {
  sl_add_val(pq_, value.first, value.second, TRANSACTIONAL);
}

bool Spraylist::Handle::try_extract_top(value_type& retval) {
  retval.first = 0;
  while (true) {
    bool success = spray_delete_min_key(pq_, &retval.first, &retval.second,
                                        data_.get()) != 0;
    if (success || retval.first == std::numeric_limits<key_type>::max())
      return success;
  }
}

}  // namespace wrapper
