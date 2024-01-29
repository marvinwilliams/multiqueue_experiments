#pragma once

#include "util.hpp"

#include "cxxopts.hpp"

#include <cstddef>
#include <mutex>
#include <optional>
#include <queue>
#include <utility>
#include <vector>

namespace wrapper::locked_pq {

template <bool Min, typename Key = unsigned long, typename T = Key>
class LockedPQ {
   public:
    using key_type = Key;
    using mapped_type = T;
    using value_type = std::pair<key_type, mapped_type>;
    using key_compare = std::conditional_t<Min, std::greater<>, std::less<>>;
    using value_compare = util::ValueCompare<value_type, key_compare, util::PairFirst>;

    struct Handle {
        friend LockedPQ;
        LockedPQ* pq_;

        Handle(LockedPQ* pq) : pq_{pq} {
        }

        void push(value_type const& value) {
            pq_->push(value);
        }

        std::optional<value_type> try_pop() {
            return pq_->try_pop();
        }
    };

   private:
    using pq_type = std::priority_queue<value_type, std::vector<value_type>, value_compare>;

    pq_type pq_{};
    std::mutex m_;

   public:
    using handle_type = util::SelfHandle<LockedPQ>;
    using settings_type = util::EmptySettings;
    LockedPQ(int /*unused*/, std::size_t capacity, settings_type const& /*unused*/) {
        std::vector<value_type> v{};
        v.reserve(capacity);
        pq_ = pq_type{value_compare{}, std::move(v)};
    }

    void push(value_type const& value) {
        std::lock_guard<std::mutex> lock{m_};
        pq_.push(value);
    }

    std::optional<value_type> try_pop() {
        std::lock_guard<std::mutex> lock{m_};
        if (pq_.empty()) {
            return std::nullopt;
        }
        auto retval = pq_.top();
        pq_.pop();
        return retval;
    }

    handle_type get_handle() {
        return handle_type{*this};
    }
};

}  // namespace wrapper::locked_pq
