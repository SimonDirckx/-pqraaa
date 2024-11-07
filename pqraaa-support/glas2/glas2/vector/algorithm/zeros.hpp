#ifndef glas2_vector_algorithm_zeros_hpp
#define glas2_vector_algorithm_zeros_hpp

#include <glas2/vector/type/static_const_vector.hpp>

namespace glas2 {

  template <typename T, typename I>
  static_const_vector<I,T,0> zeros( I n ) {
    return static_const_vector<I,T,0>( n ) ;
  }

} // namespace glas2

#endif
