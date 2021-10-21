#pragma once
#ifndef WRAPPER_KLSM_HPP_INCLUDED
#define WRAPPER_KLSM_HPP_INCLUDED

// Adapted from klsm

#include "k_lsm/k_lsm.h"

#include <limits>
#include <sstream>
#include <string>
#include <utility>

namespace wrapper {

// No known limitations

template <typename KeyType, typename T, int Relaxation>
class Klsm {
 public:
  using key_type = KeyType;
  using mapped_type = T;

  struct value_type {
    key_type key;
    mapped_type data;
  };

  struct Handle {};

  static constexpr key_type min_valid_key =
      std::numeric_limits<key_type>::min();
  static constexpr key_type max_valid_key =
      std::numeric_limits<key_type>::max();

 private:
  kpq::k_lsm<key_type, mapped_type, Relaxation> pq_;

 public:
  Klsm() = default;

  Handle get_handle() const { return Handle{}; }

  void push(Handle&, value_type value) {
    pq_.insert(value.key, value.data);
  }

  bool try_delete_min(Handle&, value_type& retval) {
    return pq_.delete_min(retval.key, retval.data);
  }

  std::string description() const {
    std::stringstream ss;
    ss << "klsm\n";
    ss << "Relaxation: " << Relaxation;
    return ss.str();
  }
};

}  // namespace wrapper

#endif
