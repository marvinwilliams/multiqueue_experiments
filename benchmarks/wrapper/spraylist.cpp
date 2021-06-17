// Adapted from klsm

#include "spraylist.hpp"

#include <cstddef>
#include <iostream>
#include <utility>

__thread unsigned long *seeds;

namespace wrapper {
thread_local thread_data_t *d;
}  // namespace wrapper
