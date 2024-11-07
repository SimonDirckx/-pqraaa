#ifndef glas2_matrix_algorithm_zeros_hpp
#define glas2_matrix_algorithm_zeros_hpp

#include <glas2/matrix/type/static_const_matrix.hpp>

namespace glas2 {

  template <typename T, typename I>
  static_const_matrix<I,T,0> zeros( I m, I n ) {
    return static_const_matrix<I,T,0>( m, n ) ;
  }

} // namespace glas2

#endif
