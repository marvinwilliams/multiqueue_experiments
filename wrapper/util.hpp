#pragma once

#include <cxxopts.hpp>

#include <optional>

namespace wrapper::util {
struct PairFirst {
    template <typename Pair>
    static constexpr auto const& get(Pair const& p) noexcept {
        return p.first;
    }
};

template <typename Value, typename KeyCompare, typename KeyOfValue>
class ValueCompare {
    [[no_unique_address]] KeyCompare comp;

   public:
    explicit ValueCompare(KeyCompare const& compare = KeyCompare{}) : comp{compare} {
    }

    constexpr bool operator()(Value const& lhs, Value const& rhs) const noexcept {
        return comp(KeyOfValue::get(lhs), KeyOfValue::get(rhs));
    }
};

template <typename PQ>
class SelfHandle {
    friend PQ;
    PQ* pq_;

    explicit SelfHandle(PQ& pq) : pq_{&pq} {
    }

   public:
    bool push(typename PQ::value_type const& value) {
        pq_->push(value);
        return true;
    }
    std::optional<typename PQ::value_type> try_pop() {
        return pq_->try_pop();
    }
};

struct EmptySettings {
    void register_cmd_options(cxxopts::Options& /*unused*/) {
    }

    [[nodiscard]] static bool validate() {
        return true;
    }

    static void write_human_readable(std::ostream& /*unused*/) {
    }

    static void write_json(std::ostream& out) {
        out << '{' << '}';
    }
};

}  // namespace wrapper::util
