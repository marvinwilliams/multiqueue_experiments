#ifndef WRAPPER_KLSM_HPP_INCLUDED
#define WRAPPER_KLSM_HPP_INCLUDED

// Adapted from klsm

#include "k_lsm/k_lsm.h"

#include <limits>
#include <sstream>
#include <string>
#include <utility>

namespace wrapper {

template <typename KeyType, typename ValueType, int Relaxation = 256>
class klsm {
 public:
  using key_type = unsigned long;
  using mapped_type = unsigned long;
  using value_type = std::pair<key_type, mapped_type>;
  struct Handle {};

  static constexpr key_type min_valid_key =
      std::numeric_limits<key_type>::min();
  static constexpr key_type max_valid_key =
      std::numeric_limits<key_type>::max();

 private:
  kpq::k_lsm<KeyType, ValueType, Relaxation> pq_;

 public:
  klsm() = default;

  Handle get_handle() const { return Handle{}; }

  void push(Handle&, value_type value) {
    pq_.insert(value.first, value.second);
  }

  bool try_delete_min(Handle&, value_type retval) {
    return pq_.delete_min(retval.first, retval.second);
  }

  std::string description() const {
    std::stringstream ss;
    ss << "klsm\n\t" << Relaxation;
    return ss.str();
  }
};

}  // namespace wrapper

#endif
