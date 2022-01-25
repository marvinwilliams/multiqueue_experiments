#pragma once
#ifndef WRAPPER_KLSM_HPP_INCLUDED
#define WRAPPER_KLSM_HPP_INCLUDED

// Adapted from klsm

#include "k_lsm/k_lsm.h"

#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <utility>

namespace wrapper {

// No known limitations

template <typename KeyType, typename T, int Relaxation>
class Klsm {
  using pq_type = kpq::k_lsm<key_type, mapped_type, Relaxation>;

 public:
  using key_type = KeyType;
  using mapped_type = T;

  using value_type = std::pair<key_type, mapped_type>;

  class Handle {
    friend Klsm;
    pq_type* pq_;

   public:
    void push(value_type const& value) {
      pq_->insert(value.first, value.second);
    }
    bool try_extract_top(value_type& retval) {
      return pq_->delete_min(retval.first, retval.second);
    }
  };

 private:
  alignas(64) std::unique_ptr<pq_type> pq_;

 public:
  Klsm(unsigned int /* num_threads */) {}

  Handle get_handle() {
    auto h = Handle{};
    h.pq_ = pq_.get();
    return h;
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
