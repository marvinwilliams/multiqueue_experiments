#pragma once

// Adapted from klsm

#include "wrapper/priority.hpp"

#include "cxxopts.hpp"

#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <ostream>
#include <utility>

namespace wrapper {

namespace detail {

class SpraylistBase {
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
    int num_threads_;

   protected:
    static constexpr key_type sentinel = std::numeric_limits<key_type>::max();

    struct ThreadDataWrapper;
    struct ThreadDataDeleter {
        void operator()(ThreadDataWrapper*);
    };

    void push(value_type const& value);
    std::optional<value_type> try_pop(ThreadDataWrapper* data);
    [[nodiscard]] std::unique_ptr<ThreadDataWrapper, ThreadDataDeleter> new_thread_data() const;

   public:
    static void add_options(cxxopts::Options& /*options*/, config_type& /*config*/) {
    }

    SpraylistBase(int num_threads, std::size_t initial_capacity, config_type const& config);

    static std::ostream& describe(std::ostream& out) {
        out << "Spraylist";
        return out;
    }
};

}  // namespace detail

template <typename Key = unsigned long, typename T = unsigned long, Priority = Priority::Min>
class Spraylist {
    static_assert(std::is_same_v<Key, unsigned long> && std::is_same_v<T, unsigned long>,
                  "Only unsigned long as Key and T are supported");
};

template <Priority P>
class Spraylist<unsigned long, unsigned long, P> : public detail::SpraylistBase {
   public:
    using detail::SpraylistBase::SpraylistBase;

    class Handle {
        friend Spraylist;
        Spraylist* pq_;
        std::unique_ptr<ThreadDataWrapper, ThreadDataDeleter> data_;

        explicit Handle(Spraylist* pq) : pq_{pq}, data_{pq_->new_thread_data()} {
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
            auto ret = pq_->try_pop(data_.get());
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
