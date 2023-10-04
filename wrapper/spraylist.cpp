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
#include <mutex>

__thread unsigned long* seeds;

namespace wrapper {

struct detail::SpraylistBase::PQWrapper : ::sl_intset_t {};

void detail::SpraylistBase::Deleter::operator()(PQWrapper* p) {
    ::sl_set_delete(p);
}

struct detail::SpraylistBase::ThreadDataWrapper : ::thread_data_t {};

void detail::SpraylistBase::ThreadDataDeleter::operator()(ThreadDataWrapper* p) {
    delete static_cast<thread_data_t*>(p);
}

detail::SpraylistBase::SpraylistBase(int num_threads, std::size_t initial_capacity,
                                                        config_type const& /*options*/)
    : num_threads_(num_threads) {
    ::ssalloc_init(num_threads_);
    *levelmax = floor_log_2(initial_capacity);
    pq_.reset(static_cast<PQWrapper*>(sl_set_new()));
}

auto detail::SpraylistBase::new_thread_data() const -> std::unique_ptr<ThreadDataWrapper, ThreadDataDeleter> {
    static std::mutex m;
    auto l = std::scoped_lock(m);
    ::ssalloc_init(num_threads_);
    seeds = ::seed_rand();
    auto thread_data = std::unique_ptr<ThreadDataWrapper, ThreadDataDeleter>(static_cast<ThreadDataWrapper*>(new thread_data_t));
    thread_data->seed = ::rand();
    thread_data->seed2 = ::rand();
    thread_data->nb_threads = num_threads_;
    return thread_data;
}

void detail::SpraylistBase::push(value_type const& value) {
    ::sl_add_val(pq_.get(), value.first, value.second, TRANSACTIONAL);
}

auto detail::SpraylistBase::try_pop(ThreadDataWrapper* data) -> std::optional<value_type> {
    value_type retval;
    int ret;
    do {
        ret = ::spray_delete_min_key(pq_.get(), &retval.first, &retval.second, data);
    } while (ret == 0 && retval.first != sentinel);
    if (retval.first == sentinel) {
        return std::nullopt;
    }
    return retval;
}

}  // namespace wrapper
