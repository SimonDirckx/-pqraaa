#ifndef glas2_matrix_algorithm_multiply_hpp
#define glas2_matrix_algorithm_multiply_hpp

#include <glas2/concept/is.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/vector/expression/binary_operation.hpp>
#include <glas2/algorithm/ops.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/matrix/expression/matrix_vector_multiply.hpp>
#include <glas2/matrix/expression/matrix_multiply.hpp>
#include <glas2/matrix/type/diagonal_matrix.hpp>
#include <type_traits>

namespace glas2 {

  template <typename M, typename V>
  typename std::enable_if< is<DenseMatrix,M>::value && is<DenseVector,V>::value
                         , matrix_vector_multiply<M,V>
                         >::type multiply( M const& m, V const& v ) {
    return matrix_vector_multiply<M,V>( m, v ) ;
  }

  template <typename M1, typename M2>
  typename std::enable_if< is<DenseMatrix,M1>::value && is<DenseMatrix,M2>::value
                         , matrix_multiply<M1,M2>
                         >::type multiply( M1 const& m1, M2 const& m2 ) {
    return matrix_multiply<M1,M2>( m1, m2 ) ;
  }

  // Diagonal matrix

  template <typename D, typename V>
  typename std::enable_if< is<DenseVector,D>::value && is<DenseVector,V>::value
                         , binary_operation<D,V,multiplies>
                         >::type multiply( diagonal_matrix<D> const& m, V const& v ) {
    return binary_operation<D,V,multiplies>( m.diagonal(), v ) ;
  }

} // namespace glas2

#endif
