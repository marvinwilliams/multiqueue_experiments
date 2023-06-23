// Adapted from klsm

#include "wrapper/spraylist.hpp"

extern "C" {
#include "spraylist_linden/include/random.h"
#include "spraylist_linden/intset.h"
#include "ssalloc.h"
}

#undef min
#undef max

#include "cxxopts.hpp"

#include <cstddef>
#include <iostream>
#include <memory>

__thread unsigned long* seeds;

namespace wrapper {

template <bool Min>
void Spraylist<unsigned long, unsigned long, Min>::sl_intset_deleter::operator()(sl_intset_t* p) {
    sl_set_delete(p);
}

template <bool Min>
void Spraylist<unsigned long, unsigned long, Min>::thread_data_deleter::operator()(thread_data_t* p) {
    delete p;
}

template <bool Min>
Spraylist<unsigned long, unsigned long, Min>::Spraylist(int num_threads, std::size_t initial_capacity,
                                                        cxxopts::ParseResult const& /*options*/)
    : num_threads_(num_threads) {
    ssalloc_init(num_threads_);
    *levelmax = floor_log_2(initial_capacity);
    pq_.reset(sl_set_new());
}

template <bool Min>
auto Spraylist<unsigned long, unsigned long, Min>::get_handle() {
    ssalloc_init(num_threads_);
    seeds = seed_rand();
    auto thread_data = std::unique_ptr<thread_data_t, thread_data_deleter>(new thread_data_t);
    thread_data->seed = rand();
    thread_data->seed2 = rand();
    thread_data->nb_threads = num_threads_;
    return Handle{std::move(thread_data), pq_.get()};
}

template <>
void Spraylist<unsigned long, unsigned long, true>::Handle::push(value_type const& value) const {
    sl_add_val(pq_, value.first, value.second, TRANSACTIONAL);
}

template <>
void Spraylist<unsigned long, unsigned long, false>::Handle::push(value_type const& value) const {
    sl_add_val(pq_, max_valid_key - value.first, value.second, TRANSACTIONAL);
}

template <bool Min>
auto Spraylist<unsigned long, unsigned long, Min>::Handle::try_pop() -> std::optional<value_type> {
    value_type retval;
    int ret;
    do {
        ret = spray_delete_min_key(pq_, &retval.first, &retval.second, data_.get());
    } while (ret == 0 && retval.first != sentinel_);
    if (retval.first == sentinel_) {
        return std::nullopt;
    }
    if (!Min) {
        retval.first = max_valid_key - retval.first;
    }
    return retval;
}

template class Spraylist<unsigned long, unsigned long, true>;
template class Spraylist<unsigned long, unsigned long, false>;

}  // namespace wrapper
