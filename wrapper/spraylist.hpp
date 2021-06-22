#ifndef WRAPPER_SPRAYLIST_HPP_INCLUDED
#define WRAPPER_SPRAYLIST_HPP_INCLUDED

// Adapted from klsm

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>

struct sl_intset;

namespace wrapper {

class spraylist {
  using pq_t = sl_intset;

  pq_t* pq_;

 public:
  struct Handle {};
  spraylist(size_t const num_threads);
  virtual ~spraylist();

  constexpr Handle get_handle(unsigned int) { return Handle{}; }

  void init_thread(size_t const num_threads);

  /* typedef unsigned long slkey_t; */
  /* typedef unsigned long val_t; */
  /* #define KEY_MIN                         0 */
  /* #define KEY_MAX                         UINT32_MAX */
  void push(Handle, std::pair<unsigned long, unsigned long> const& value);
  bool extract_top(Handle, std::pair<unsigned long, unsigned long>& retval);

  static std::string description() { return "spraylist"; }
};

}  // namespace wrapper

#endif
