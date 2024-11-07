#ifndef glas2_matrix_algorithm_ones_hpp
#define glas2_matrix_algorithm_ones_hpp

#include <glas2/matrix/type/static_const_matrix.hpp>

namespace glas2 {

  template <typename T, typename I>
  static_const_matrix<I,T,1> ones( I m, I n ) {
    return static_const_matrix<I,T,1>( m, n ) ;
  }

} // namespace glas2

#endif
