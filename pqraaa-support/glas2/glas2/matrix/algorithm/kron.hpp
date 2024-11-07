//  (C) Copyright Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_algorithm_kron_hpp
#define glas2_matrix_algorithm_kron_hpp

#include <glas2/concept/is.hpp>
#include <glas2/matrix/expression/matrix_kron_expression.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <type_traits>

namespace glas2 {

  template <typename V, typename W>
  typename std::enable_if< is<DenseMatrix,V>::value && is<DenseMatrix,W>::value, matrix_kron_expression< V, W > >::type kron( V const& v, W const& w ) {
    return matrix_kron_expression< V, W >( v, w ) ;
  }

} // namespace glas2

#endif
