#pragma once
#ifndef WRAPPER_CAPQ_HPP_INCLUDED
#define WRAPPER_CAPQ_HPP_INCLUDED

// Adapted from klsm

#include <ios>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <utility>

namespace wrapper {

// GC has MAX_THREADS = 128
// (unsigned long) -1 is signaling empty

template <bool remove_min_relax = true, bool put_relax = true,
          bool catree_adapt = true>
class Capq {
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
  static constexpr key_type max_valid_key =
      std::numeric_limits<key_type>::max() - 1;

 private:
  static constexpr key_type empty_key = std::numeric_limits<key_type>::max();

  struct wrapper_type;

  std::unique_ptr<wrapper_type> pq_;

 public:
  Capq();
  ~Capq();

  Handle get_handle() const { return Handle{}; }

  void push(Handle&, value_type value);
  bool try_delete_min(Handle&, value_type& retval);

  std::string description() const {
    std::stringstream ss;
    ss << "capq\n";
    ss << std::boolalpha;
    ss << "Remove min relax: " << remove_min_relax << '\n';
    ss << "Put relax: " << put_relax << '\n';
    ss << "Catree adapt: " << catree_adapt << '\n';
    return ss.str();
  }
};

}  // namespace wrapper

#endif
