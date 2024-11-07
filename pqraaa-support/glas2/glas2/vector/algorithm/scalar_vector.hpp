//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_algorithm_scalar_vector_hpp
#define glas2_vector_algorithm_scalar_vector_hpp

#include <glas2/vector/expression/scalar_vector_expression.hpp>
#include <glas2/algorithm/ops.hpp>
#include <glas2/concept/multiplies.hpp>
#include <glas2/concept/plus.hpp>
#include <glas2/concept/minus.hpp>
#include <type_traits>

namespace glas2 {

/*  template <typename S, typename V>
  typename std::enable_if< is<Vector,V>::value, scalar_vector_expression<S,V,glas2::multiplies> >::type operator*( S const& s, V const& v ) {
    return scalar_vector_expression<S,V,glas2::multiplies>( s, v ) ;
  }

  template <typename S, typename V>
  typename std::enable_if< is<Vector,V>::value, scalar_vector_expression<S,V,glas2::plus> >::type operator+( S const& s, V const& v ) {
    return scalar_vector_expression<S,V,glas2::plus>( s, v ) ;
  }

  template <typename S, typename V>
  typename std::enable_if< is<Vector,V>::value, scalar_vector_expression<S,V,glas2::minus> >::type operator-( S const& s, V const& v ) {
    return scalar_vector_expression<S,V,glas2::minus>( s, v ) ;
  }*/

} // namespace glas2

#endif
