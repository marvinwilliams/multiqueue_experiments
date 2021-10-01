#ifndef WRAPPER_CAPQ_HPP_INCLUDED
#define WRAPPER_CAPQ_HPP_INCLUDED

// Adapted from klsm

extern "C" {
#include "capq/capq.h"
#include "capq/gc/gc.h"
}

#include <cstddef>
#include <cstdint>
#include <limits>
#include <sstream>
#include <string>
#include <utility>

namespace wrapper {

template <bool remove_min_relax = true, bool put_relax = true,
          bool catree_adapt = true>
class capq {
 public:
  using key_type = unsigned long;
  using mapped_type = unsigned long;
  using value_type = std::pair<key_type, mapped_type>;
  struct Handle {};

  static constexpr key_type min_valid_key =
      std::numeric_limits<key_type>::min();
  static constexpr key_type max_valid_key =
      std::numeric_limits<key_type>::max() - 1;

 private:
  struct capq_wrapper_type {
    char pad1[64 - sizeof(CAPQ*)];
    CAPQ* pq;
    char pad2[64];
  };

  static constexpr key_type empty_key = std::numeric_limits<key_type>::max();

  capq_wrapper_type* pq_;

 public:
  capq() {
    _init_gc_subsystem();
    pq_ = new capq_wrapper_type;
    pq_->pq = capq_new();
  }

  ~capq() {
    capq_delete(pq_->pq);
    delete pq_;
    _destroy_gc_subsystem();
  }

  Handle get_handle() const { return Handle{}; }

  void push(Handle&, value_type value) {
    capq_put_param(pq_->pq, value.first, value.second, catree_adapt);
  }
  bool try_delete_min(Handle&, value_type& retval) {
    retval.second = capq_remove_min_param(
        pq_->pq, &retval.first, remove_min_relax, put_relax, catree_adapt);
    return retval.first != empty_key;
  }

  std::string description() const {
    std::stringstream ss;
    ss << "capq\n";
    ss << "Remove min relax: " << std::boolalpha << remove_min_relax << '\n';
    ss << "Put relax: " << std::boolalpha << put_relax << '\n';
    ss << "Catree adapt: " << std::boolalpha << catree_adapt << '\n';
    return ss.str();
  }
};

template class capq<true, true, true>;
template class capq<true, false, true>;
template class capq<false, true, true>;
template class capq<false, false, true>;


}  // namespace wrapper

#endif
