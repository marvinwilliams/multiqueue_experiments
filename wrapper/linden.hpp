#pragma once
#ifndef WRAPPER_LINDEN_HPP_INCLUDED
#define WRAPPER_LINDEN_HPP_INCLUDED

// Adapted from klsm

#include <limits>
#include <memory>
#include <string>
#include <utility>

namespace wrapper {

class Linden {
 public:
  using key_type = unsigned long;
  using mapped_type = unsigned long;

  struct value_type {
    key_type key;
    mapped_type data;
  };

  struct Handle {};

  static constexpr key_type min_valid_key =
      std::numeric_limits<key_type>::min();  // The queue itself only supports
                                             // keys >= 1, so one is added on
                                             // each insert
  static constexpr key_type max_valid_key =
      std::numeric_limits<key_type>::max() - 2;

 private:
  static constexpr key_type empty_key =
      std::numeric_limits<key_type>::max() - 1;

  struct wrapper_type;

  alignas(64) std::unique_ptr<wrapper_type> pq_;

 public:
  Linden();
  ~Linden();

  Handle get_handle() { return Handle{}; }

  void push(Handle&, value_type value);

  bool try_delete_min(Handle&, value_type& retval);

  std::string description() const { return "linden"; }
};

}  // namespace wrapper

#endif
