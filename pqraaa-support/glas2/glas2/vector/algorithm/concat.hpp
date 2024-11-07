//  (C) Copyright Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_algorithm_concat_hpp
#define glas2_vector_algorithm_concat_hpp

#include <glas2/concept/is.hpp>
#include <glas2/vector/expression/vector_concat_expression.hpp>
#include <glas2/vector/expression/scalar_concat_expression.hpp>
#include <glas2/vector/expression/scalar_vector_concat_expression.hpp>
#include <glas2/vector/expression/vector_scalar_concat_expression.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <type_traits>

namespace glas2 {

  template <typename V, typename W>
  typename std::enable_if< is<DenseVector,V>::value && is<DenseVector,W>::value, vector_concat_expression< V, W > >::type concat( V const& v, W const& w ) {
    return vector_concat_expression< V, W >( v, w ) ;
  }

  template <typename V, typename W>
  typename std::enable_if< is<Scalar,V>::value && is<DenseVector,W>::value, scalar_vector_concat_expression< V, W > >::type concat( V const& v, W const& w ) {
    return scalar_vector_concat_expression< V, W >( v, w ) ;
  }

  template <typename V, typename W>
  typename std::enable_if< is<DenseVector,V>::value && is<Scalar,W>::value, vector_scalar_concat_expression< V, W > >::type concat( V const& v, W const& w ) {
    return vector_scalar_concat_expression< V, W >( v, w ) ;
  }

  template <typename V, typename W>
  typename std::enable_if< is<Scalar,V>::value && is<Scalar,W>::value, scalar_concat_expression< V, W > >::type concat( V const& v, W const& w ) {
    return scalar_concat_expression< V, W >( v, w ) ;
  }

} // namespace glas2

#endif
