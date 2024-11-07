//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_algorithm_scalar_matrix_hpp
#define glas2_matrix_algorithm_scalar_matrix_hpp

#include <glas2/matrix/expression/scalar_matrix_expression.hpp>
#include <glas2/concept/multiplies.hpp>
#include <glas2/concept/plus.hpp>
#include <glas2/concept/minus.hpp>
#include <type_traits>

namespace glas2 {

  template <typename S, typename M>
  typename std::enable_if< is<Matrix,M>::value, scalar_matrix_expression<S,M,glas2::multiplies> >::type operator*( S const& s, M const& m ) {
    return scalar_matrix_expression<S,M,glas2::multiplies>( s, m ) ;
  }

  template <typename S, typename M>
  typename std::enable_if< is<Matrix,M>::value, scalar_matrix_expression<S,M,glas2::plus> >::type operator+( S const& s, M const& m ) {
    return scalar_matrix_expression<S,M,glas2::plus>( s, m ) ;
  }

  template <typename S, typename M>
  typename std::enable_if< is<Matrix,M>::value, scalar_matrix_expression<S,M,glas2::minus> >::type operator-( S const& s, M const& m ) {
    return scalar_matrix_expression<S,M,glas2::minus>( s, m ) ;
  }

} // namespace glas2

#endif
