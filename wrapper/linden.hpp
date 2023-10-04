#pragma once

// Adapted from klsm

#include "wrapper/priority.hpp"

#include "cxxopts.hpp"

#include <limits>
#include <memory>
#include <optional>
#include <ostream>
#include <utility>

namespace wrapper {

namespace detail {

// The queue itself only supports
// keys >= 1, so one is added on each insert
class LindenBase {
   public:
    using key_type = unsigned long;
    using mapped_type = unsigned long;
    using value_type = std::pair<key_type, mapped_type>;
    struct config_type {};

   private:
    struct PQWrapper;
    struct Deleter {
        void operator()(PQWrapper*);
    };

    alignas(64) std::unique_ptr<PQWrapper, Deleter> pq_;

   protected:
    static constexpr key_type sentinel = std::numeric_limits<key_type>::max();

    void push(value_type const& value);
    std::optional<value_type> try_pop();

   public:
    static void add_options(cxxopts::Options& /*options*/, config_type& /*config*/) {
    }

    LindenBase(int num_threads, std::size_t initial_capacity, config_type const& options);

    static std::ostream& describe(std::ostream& out) {
        out << "Linden";
        return out;
    }
};

}  // namespace detail

template <typename Key = unsigned long, typename T = unsigned long, Priority = Priority::Min>
class Linden {
    static_assert(std::is_same_v<Key, unsigned long> && std::is_same_v<T, unsigned long>,
                  "Only unsigned long as Key and T are supported");
};

template <Priority P>
class Linden<unsigned long, unsigned long, P> : public detail::LindenBase {
   public:
    using detail::LindenBase::LindenBase;

    class Handle {
        friend Linden;
        Linden* pq_;
        explicit Handle(Linden* pq) : pq_{pq} {
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
