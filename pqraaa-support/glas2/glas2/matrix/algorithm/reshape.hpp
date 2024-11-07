#ifndef glas2_matrix_algorithm_reshape_hpp
#define glas2_matrix_algorithm_reshape_hpp

#include <glas2/concept/is.hpp>
#include <glas2/vector/concept/contiguous_dense_vector.hpp>
#include <glas2/matrix/type/contiguous_matrix.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename M, typename I1, typename I2>
  typename std::enable_if< is<ContiguousDenseMatrix,M>::value
                         , contiguous_matrix<typename M::value_type,typename M::size_type,typename M::orientation>
                         >::type reshape( M m, I1 rows, I2 columns ) {
    assert( rows*columns==m.num_rows()*m.num_columns() ) ;
    return contiguous_matrix<typename M::value_type,typename M::size_type,typename M::orientation>( m.ptr(), rows, columns ) ;
  }

  template <typename V, typename I1, typename I2, typename O>
  typename std::enable_if< is<ContiguousDenseVector,V>::value
                         , contiguous_matrix<typename V::value_type,typename V::size_type,O>
                         >::type reshape( V v, I1 rows, I2 columns, O ) {
    assert( rows*columns==v.size() ) ;
    return contiguous_matrix<typename V::value_type,typename V::size_type,O>( v.ptr(), rows, columns ) ;
  }

} // namespace glas2

#endif
