#pragma once
#ifndef WRAPPER_SPRAYLIST_HPP_INCLUDED
#define WRAPPER_SPRAYLIST_HPP_INCLUDED

// Adapted from klsm

#include <cstdint>
#include <limits>
#include <memory>
#include <ostream>
#include <utility>

using sl_intset_t = struct sl_intset;
using thread_data_t = struct thread_data;

namespace wrapper {

class Spraylist {
    struct sl_intset_deleter {
        void operator()(sl_intset_t*);
    };
    struct thread_data_deleter {
        void operator()(thread_data_t*);
    };

   public:
    using key_type = unsigned long;
    using mapped_type = unsigned long;
    using value_type = std::pair<key_type, mapped_type>;

    struct Handle {
        std::unique_ptr<thread_data_t, thread_data_deleter> data_;
        sl_intset_t* pq_;

        void push(value_type const& value) const;
        bool try_pop(value_type& retval) const;
    };

    using handle_type = Handle;

    static constexpr key_type min_valid_key = std::numeric_limits<key_type>::min();
    // Use INT_MAX_32
    static constexpr key_type max_valid_key = std::numeric_limits<std::uint32_t>::max() - 1UL;
    static_assert(std::numeric_limits<unsigned long>::max() == -1UL);

   private:
    static constexpr key_type sentinel_ = std::numeric_limits<key_type>::max();

    alignas(64) std::unique_ptr<sl_intset_t, sl_intset_deleter> pq_;
    int num_threads_;
    std::unique_ptr<thread_data_t, thread_data_deleter> thread_data_;

   public:
    explicit Spraylist(int num_threads);

    Handle get_handle();

    void push(value_type const& value) const;
    bool try_pop(value_type& retval) const;

    static std::ostream& describe(std::ostream& out) {
        out << "spraylist\n";
        return out;
    }
};

}  // namespace wrapper

#endif
