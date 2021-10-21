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

__thread unsigned long *seeds;

namespace wrapper {

struct Spraylist::wrapper_type {
  sl_intset_t *pq;
  ~wrapper_type() { sl_set_delete(pq); }
};

static thread_local std::unique_ptr<thread_data_t> d;

Spraylist::Spraylist() : pq_(new wrapper_type) {
  static constexpr unsigned int INITIAL_SIZE = 1000000;
  *levelmax = floor_log_2(INITIAL_SIZE);
  init_thread(1);
  pq_->pq = sl_set_new();
}

Spraylist::~Spraylist() {}

void Spraylist::init_thread(size_t num_threads) {
  if (!d) {
    ssalloc_init(num_threads);
    seeds = seed_rand();

    d.reset(new thread_data_t);
    d->seed = rand();
    d->seed2 = rand();
  }
  d->nb_threads = num_threads;
}

void Spraylist::push(Handle &, value_type value) {
  sl_add_val(pq_->pq, value.key, value.data, TRANSACTIONAL);
}

bool Spraylist::try_delete_min(Handle &, value_type &retval) {
  do {
  } while (spray_delete_min_key(pq_->pq, &retval.key, &retval.data,
                                d.get()) == 0 &&
           retval.key != empty_key);
  return retval.key != empty_key;
}

}  // namespace wrapper
