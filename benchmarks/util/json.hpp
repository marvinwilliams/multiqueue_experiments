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
        // std::filesystem::path quotes itself, arithmetic types stream as-is.
        // Stream format flags (e.g. std::fixed) are the caller's business.
        out << value;
    }
}

}  // namespace detail

// Writes a JSON object and keeps track of the separators, which is the only
// part of emitting these result files that is easy to get wrong. The closing
// brace is written by the destructor, so an Object must go out of scope before
// anything else is written to the same stream.
class Object {
    std::ostream* out_;
    bool empty_{true};

    std::ostream& key(std::string_view name) {
        if (!empty_) {
            *out_ << ',';
        }
        empty_ = false;
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

    //! Nested object; the callable receives the writer for it
    template <typename F>
    Object& object(std::string_view name, F&& write_body) {
        Object nested{key(name)};
        std::forward<F>(write_body)(nested);
        return *this;
    }

    //! Array of elements; the callable writes a single element
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

    //! Array of scalars
    template <typename Range>
    Object& array(std::string_view name, Range const& range) {
        return array(name, std::begin(range), std::end(range),
                     [](std::ostream& out, auto const& e) { detail::write_value(out, e); });
    }

    //! Escape hatch for values written by code that only knows about streams
    template <typename F>
    Object& raw(std::string_view name, F&& write_raw) {
        std::forward<F>(write_raw)(key(name));
        return *this;
    }
};

}  // namespace json
