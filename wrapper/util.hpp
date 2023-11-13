#pragma once

#include <type_traits>
#include <utility>

namespace wrapper::util {
template <typename Key, typename Value>
struct KeyOfValue {
    static_assert(std::is_same_v<Key, Value>, "KeyOfValue not specialized for this value type");
};

template <typename Key>
struct KeyOfValue<Key, Key> {
    static constexpr Key const& get(Key const& key) noexcept {
        return key;
    }
};

template <typename Key, typename T>
struct KeyOfValue<Key, std::pair<Key, T>> {
    static constexpr Key const& get(std::pair<Key, T> const& p) noexcept {
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

}  // namespace wrapper::util
