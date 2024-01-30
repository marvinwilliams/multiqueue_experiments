#pragma once

#include "wrapper/priority.hpp"

#include "util.hpp"

#include <cxxopts.hpp>

#include <tbb/concurrent_queue.h>

#include <optional>
#include <ostream>
#include <utility>

namespace wrapper::tbb_fifo {

// Min is ignored since this is a FIFO
template <bool /*Min*/, typename Key = unsigned long, typename T = Key>
class TBBFIFO {
   public:
    using key_type = Key;
    using mapped_type = T;
    using value_type = std::pair<key_type, mapped_type>;

   private:
    using pq_type = tbb::concurrent_queue<value_type>;

    pq_type pq_;

   public:
    using handle_type = util::SelfHandle<TBBFIFO>;
    using settings_type = util::EmptySettings;
    TBBFIFO(int /*num_threads*/, std::size_t /*initial_capacity*/, settings_type const& /*options*/) {
    }

    void push(value_type const& value) {
        pq_.push(value);
    }
    std::optional<value_type> try_pop() {
        value_type retval;
        if (!pq_.try_pop(retval)) {
            return std::nullopt;
        }
        return retval;
    }

    static void write_human_readable(std::ostream& out) {
        out << "TBB FIFO";
    }

    handle_type get_handle() {
        return handle_type{*this};
    }
};

}  // namespace wrapper
