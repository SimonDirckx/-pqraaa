//  (C) Copyright Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_algorithm_row_kron_hpp
#define glas2_matrix_algorithm_row_kron_hpp

#include <glas2/concept/is.hpp>
#include <glas2/matrix/expression/vector_matrix_row_kron_expression.hpp>
#include <glas2/matrix/expression/matrix_vector_row_kron_expression.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/matrix/concept/dense_vector.hpp>
#include <type_traits>

namespace glas2 {

  template <typename V, typename W>
  typename std::enable_if< is<DenseVector,V>::value && is<DenseMatrix,W>::value, vector_matrix_row_kron_expression< V, W > >::type row_kron( V const& v, W const& w ) {
    return vector_matrix_row_kron_expression< V, W >( v, w ) ;
  }

  template <typename V, typename W>
  typename std::enable_if< is<DenseMatrix,V>::value && is<DenseVector,W>::value, matrix_vector_row_kron_expression< V, W > >::type row_kron( V const& v, W const& w ) {
    return matrix_vector_row_kron_expression< V, W >( v, w ) ;
  }

} // namespace glas2

#endif
