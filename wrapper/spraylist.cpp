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

namespace wrapper {

void Spraylist::sl_intset_deleter::operator()(sl_intset_t* p) {
    sl_set_delete(p);
}
void Spraylist::thread_data_deleter::operator()(thread_data_t* p) {
    delete p;
}

Spraylist::Spraylist(int num_threads, std::size_t initial_capacity) : num_threads_(num_threads) {
    ssalloc_init(num_threads_);
    *levelmax = floor_log_2(initial_capacity);
    pq_.reset(sl_set_new());
}

Spraylist::Handle Spraylist::get_handle(int /*unused*/) {
    ssalloc_init(num_threads_);
    seeds = seed_rand();
    auto thread_data = std::unique_ptr<thread_data_t, thread_data_deleter>(new thread_data_t);
    thread_data->seed = rand();
    thread_data->seed2 = rand();
    thread_data->nb_threads = num_threads_;
    return Handle{std::move(thread_data), pq_.get()};
}

void Spraylist::Handle::push(value_type const& value) const {
    sl_add_val(pq_, value.first, value.second, TRANSACTIONAL);
}

bool Spraylist::Handle::try_pop(value_type& retval) const {
    int ret;
    do {
        ret = spray_delete_min_key(pq_, &retval.first, &retval.second, data_.get());
    } while (ret == 0 && retval.first != std::numeric_limits<key_type>::max());
    return retval.first != std::numeric_limits<key_type>::max();
}

}  // namespace wrapper
