#pragma once
#include <cstdint>
#ifndef WRAPPER_SPRAYLIST_HPP_INCLUDED
#define WRAPPER_SPRAYLIST_HPP_INCLUDED

// Adapted from klsm

#include <limits>
#include <memory>
#include <string>
#include <utility>

namespace wrapper {

class Spraylist {
 public:
  using key_type = unsigned long;
  using mapped_type = unsigned long;
  struct value_type {
    key_type key;
    mapped_type data;
  };

  struct Handle {};

  static constexpr key_type min_valid_key =
      std::numeric_limits<key_type>::min();
  static constexpr key_type max_valid_key = std::numeric_limits<std::uint32_t>::max() - 1;

 private:
  static constexpr key_type empty_key = std::numeric_limits<key_type>::max();

  struct wrapper_type;

  alignas(64) std::unique_ptr<wrapper_type> pq_;

 public:
  Spraylist();
  ~Spraylist();

  Handle get_handle() { return Handle{}; }

  void init_thread(size_t const num_threads);

  void push(Handle&, value_type value);
  bool try_delete_min(Handle&, value_type& retval);

  std::string description() const { return "spraylist"; }
};

}  // namespace wrapper

#endif
