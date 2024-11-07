#ifndef cork_utility_is_infinite_hpp
#define cork_utility_is_infinite_hpp

#include <limits>
#include <cmath>
#include <complex>

namespace CORK {

  template <typename V>
  bool is_infinite( V const v ) {
    return std::numeric_limits<decltype(std::abs(V()))>::infinity()==v ;
  }

} // namespace CORK

#endif
