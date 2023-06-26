#pragma once

#include <optional>

struct OperationCounters {
    long long pop = 0;
    long long push = 0;
    long long failed_pop = 0;
};

template <typename PriorityQueue>
class CountingHandle : public PriorityQueue::handle_type {
    using base_type = typename PriorityQueue::handle_type;
    OperationCounters counters_;

   public:
    explicit CountingHandle(PriorityQueue& pq) : base_type(pq.get_handle()) {
    }

    void push(typename PriorityQueue::key_type const& key) {
        base_type::push({key, key});
        ++counters_.push;
    }

    std::optional<typename PriorityQueue::value_type> try_pop() {
        auto ret = base_type::try_pop();
        if (!ret) {
            ++counters_.failed_pop;
            return std::nullopt;
        }
        ++counters_.pop;
        return *ret;
    }

    auto const& get_operation_counts() const {
        return counters_;
    }
};
