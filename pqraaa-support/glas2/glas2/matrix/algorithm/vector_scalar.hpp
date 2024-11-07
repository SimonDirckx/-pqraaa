//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_algorithm_vector_scalar_hpp
#define glas2_vector_algorithm_vector_scalar_hpp

#include <glas2/vector/expression/vector_scalar_expression.hpp>
#include <glas2/concept/divides.hpp>
#include <glas2/concept/multiplies.hpp>
#include <glas2/concept/plus.hpp>
#include <glas2/concept/minus.hpp>
#include <type_traits>

namespace glas2 {

  template <typename V, typename S>
  typename std::enable_if< is<Vector,V>::value, vector_scalar_expression<V,S,glas2::divides> >::type operator/( V const& v, S const& s ) {
    return vector_scalar_expression<V,S,glas2::divides>( v, s ) ;
  }

  template <typename V, typename S>
  typename std::enable_if< is<Vector,V>::value, vector_scalar_expression<V,S,glas2::multiplies> >::type operator*( V const& v, S const& s ) {
    return vector_scalar_expression<V,S,glas2::multiplies>( v, s ) ;
  }

  template <typename V, typename S>
  typename std::enable_if< is<Vector,V>::value, vector_scalar_expression<V,S,glas2::plus> >::type operator+( V const& v, S const& s ) {
    return vector_scalar_expression<V,S,glas2::plus>( v, s ) ;
  }

  template <typename V, typename S>
  typename std::enable_if< is<Vector,V>::value, vector_scalar_expression<V,S,glas2::minus> >::type operator-( V const& v, S const& s ) {
    return vector_scalar_expression<V,S,glas2::minus>( v, s ) ;
  }

} // namespace glas2

#endif
