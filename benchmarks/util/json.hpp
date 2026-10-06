#pragma once

#include <iomanip>
#include <iterator>
#include <ostream>
#include <string_view>
#include <type_traits>
#include <utility>

namespace json {

namespace detail {

template <typename T>
void write_value(std::ostream& out, T const& value) {
    if constexpr (std::is_convertible_v<T const&, std::string_view>) {
        out << std::quoted(std::string_view{value});
    } else {
        out << value;
    }
}

}  // namespace detail

class Object {
    std::ostream* out_;
    bool first_{true};

    std::ostream& key(std::string_view name) {
        if (!first_) {
            *out_ << ',';
        }
        first_ = false;
        *out_ << std::quoted(name) << ':';
        return *out_;
    }

   public:
    explicit Object(std::ostream& out) : out_{&out} {
        *out_ << '{';
    }

    Object(Object const&) = delete;
    Object(Object&&) = delete;
    Object& operator=(Object const&) = delete;
    Object& operator=(Object&&) = delete;

    ~Object() {
        *out_ << '}';
    }

    template <typename T>
    Object& entry(std::string_view name, T const& value) {
        detail::write_value(key(name), value);
        return *this;
    }

    template <typename F>
    Object& object(std::string_view name, F&& write_inner) {
        Object inner{key(name)};
        std::forward<F>(write_inner)(inner);
        return *this;
    }

    template <typename It, typename F>
    Object& array(std::string_view name, It first, It last, F write_element) {
        auto& out = key(name);
        out << '[';
        for (auto it = first; it != last; ++it) {
            if (it != first) {
                out << ',';
            }
            write_element(out, *it);
        }
        out << ']';
        return *this;
    }

    template <typename Range>
    Object& array(std::string_view name, Range const& range) {
        return array(name, std::begin(range), std::end(range),
                     [](std::ostream& out, auto const& e) { detail::write_value(out, e); });
    }

    template <typename F>
    Object& raw(std::string_view name, F&& write_raw) {
        std::forward<F>(write_raw)(key(name));
        return *this;
    }
};

}  // namespace json

