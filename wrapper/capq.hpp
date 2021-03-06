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

typedef struct fpasl_catree_set CAPQ;
struct CAPQ_deleter {
  void operator()(CAPQ*);
};

namespace wrapper {

// GC has MAX_THREADS = 128
// (unsigned long) -1 is signaling empty

template <bool remove_min_relax = true, bool put_relax = true,
          bool catree_adapt = true>
class Capq {
 public:
  using key_type = unsigned long;
  using mapped_type = unsigned long;
  using value_type = std::pair<key_type, mapped_type>;

  struct Handle {
    CAPQ* pq_;

    void push(value_type const& value);
    bool try_pop(value_type& retval);
  };

  static constexpr key_type min_valid_key =
      std::numeric_limits<key_type>::min();
  static constexpr key_type max_valid_key =
      std::numeric_limits<key_type>::max() - 1;

 private:
  static constexpr key_type sentinel_ = std::numeric_limits<key_type>::max();

  alignas(64) std::unique_ptr<CAPQ, CAPQ_deleter> pq_;

 public:
  Capq(unsigned int /* num_threads */);

  Handle get_handle();

  void push(value_type const& value);
  bool try_pop(value_type& retval);

  std::string description() const {
    std::stringstream ss;
    ss << "capq\n";
    ss << std::boolalpha;
    ss << "Remove min relax: " << remove_min_relax << '\n';
    ss << "Put relax: " << put_relax << '\n';
    ss << "Catree adapt: " << catree_adapt;
    return ss.str();
  }
};

}  // namespace wrapper

#endif
