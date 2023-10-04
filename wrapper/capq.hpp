#pragma once

// Adapted from klsm

#include "wrapper/priority.hpp"

#include "cxxopts.hpp"

#include <ios>
#include <limits>
#include <memory>
#include <optional>
#include <ostream>
#include <utility>

using capq_type = struct fpasl_catree_set;

namespace wrapper {

namespace detail {
// GC has MAX_THREADS = 128
// (unsigned long) -1 is signaling empty

class CAPQBase {
   public:
    using key_type = unsigned long;
    using mapped_type = unsigned long;
    using value_type = std::pair<key_type, mapped_type>;
    struct config_type {};

   private:
    struct Deleter {
        void operator()(capq_type*);
    };
    alignas(64) std::unique_ptr<capq_type, Deleter> pq_;

   protected:
    static constexpr key_type sentinel = std::numeric_limits<key_type>::max();
    void push(value_type const& value);
    std::optional<value_type> try_pop();

   public:
    static void add_options(cxxopts::Options& /*options*/, config_type& /*config*/) {
    }

    CAPQBase(int num_threads, std::size_t initial_capacity, config_type const& config);

    static std::ostream& describe(std::ostream& out) {
        out << "CA-PQ";
        return out;
    }
};

}  // namespace detail

template <typename Key = unsigned long, typename T = unsigned long, Priority = Priority::Min>
class CAPQ {
    static_assert(std::is_same_v<Key, unsigned long> && std::is_same_v<T, unsigned long>,
                  "Only unsigned long as Key and T are supported");
};

template <Priority P>
class CAPQ<unsigned long, unsigned long, P> : public detail::CAPQBase {
   public:
    using detail::CAPQBase::CAPQBase;

    class Handle {
        friend CAPQ;
        CAPQ* pq_;

        explicit Handle(CAPQ* pq) : pq_{pq} {
        }

       public:
        void push(value_type const& value) {
            if constexpr (P == Priority::Max) {
                pq_->push({sentinel - value.first - 1, value.second});
            } else {
                pq_->push(value);
            }
        }

        std::optional<value_type> try_pop() {
            auto ret = pq_->try_pop();
            if constexpr (P == Priority::Max) {
                if (ret) {
                    ret->first = sentinel - ret->first - 1;
                }
            }
            return ret;
        };
    };

    using handle_type = Handle;

    Handle get_handle() {
        return Handle{this};
    }
};

}  // namespace wrapper
