#ifndef WRAPPER_LINDEN_HPP_INCLUDED
#define WRAPPER_LINDEN_HPP_INCLUDED

// Adapted from klsm

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>

namespace wrapper {

struct linden_pq_wrapper;

class linden {
  linden_pq_wrapper* pq_;

 public:
  struct Handle {};
  static constexpr int DEFAULT_OFFSET = 32;

  linden(unsigned int num_threads = 0, int const max_offset = DEFAULT_OFFSET);

  constexpr Handle get_handle(unsigned int) { return Handle{}; }

  ~linden();

  /* #define SENTINEL_KEYMIN (0UL)  /1* Key value of first dummy node. *1/ */
  /* #define SENTINEL_KEYMAX (~1UL) /1* Key value of last dummy node.  *1/ */
  void push(Handle, std::pair<unsigned long, unsigned long> const& value);

  bool extract_top(Handle, std::pair<unsigned long, unsigned long>& retval);

  static std::string description() { return "linden"; }
};

}  // namespace wrapper

#endif
