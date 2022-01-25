#pragma once
#ifndef WRAPPER_SPRAYLIST_HPP_INCLUDED
#define WRAPPER_SPRAYLIST_HPP_INCLUDED

// Adapted from klsm

extern "C" {
#include "spraylist_linden/intset.h"
}

#undef min
#undef max
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <utility>

namespace wrapper {

class Spraylist {
 public:
  using key_type = unsigned long;
  using mapped_type = unsigned long;
  using value_type = std::pair<key_type, mapped_type>;

  class Handle {
    friend Spraylist;
    std::unique_ptr<thread_data_t> data_;
    sl_intset_t* pq_;

   public:
    void push(value_type const& value);
    bool try_extract_top(value_type& retval);
  };

  static constexpr key_type min_key = std::numeric_limits<key_type>::min();
  // Use INT_MAX_32
  static constexpr key_type max_key =
      std::numeric_limits<std::uint32_t>::max() - 1;
  static_assert(std::numeric_limits<unsigned long>::max() == -1, "");

 private:
  static constexpr key_type sentinel_ = std::numeric_limits<key_type>::max();

  alignas(64) std::unique_ptr<sl_intset_t, void (*)(sl_intset_t*)> pq_;
  unsigned int num_threads_;

 public:
  explicit Spraylist(unsigned int num_threads);

  Handle get_handle();

  std::string description() const { return "spraylist"; }
};

}  // namespace wrapper

#endif
