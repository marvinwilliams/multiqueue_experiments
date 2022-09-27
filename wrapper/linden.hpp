#pragma once
#ifndef WRAPPER_LINDEN_HPP_INCLUDED
#define WRAPPER_LINDEN_HPP_INCLUDED

// Adapted from klsm

#include <limits>
#include <memory>
#include <ostream>
#include <utility>

struct pq_t;

namespace wrapper {

class Linden {
    struct pq_deleter {
        void operator()(pq_t*);
    };

   public:
    using key_type = unsigned long;
    using mapped_type = unsigned long;

    using value_type = std::pair<key_type, mapped_type>;

    struct Handle {
        pq_t* pq_;

        void push(value_type const& value) const;
        bool try_pop(value_type& retval) const;
    };

    using handle_type = Handle;

    // The queue itself only supports
    // keys >= 1, so one is added on each insert
    static constexpr key_type min_valid_key = std::numeric_limits<key_type>::min();
    static constexpr key_type max_valid_key = std::numeric_limits<key_type>::max() - 3;

   private:
    static constexpr key_type sentinel_ = std::numeric_limits<key_type>::max();

    struct wrapper_type;

    alignas(64) std::unique_ptr<pq_t, pq_deleter> pq_;

   public:
    Linden(unsigned int num_threads);

    Handle get_handle();

    void push(value_type const& value) const;
    bool try_pop(value_type& retval) const;

    static std::ostream& describe(std::ostream& out) {
        out << "linden\n";
        return out;
    }
};

}  // namespace wrapper

#endif
