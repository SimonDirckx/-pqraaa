//  (C) Copyright Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_algorithm_kron_hpp
#define glas2_vector_algorithm_kron_hpp

#include <glas2/concept/is.hpp>
#include <glas2/vector/expression/vector_kron_expression.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <type_traits>

namespace glas2 {

  template <typename V, typename W>
  typename std::enable_if< is<DenseVector,V>::value && is<DenseVector,W>::value, vector_kron_expression< V, W > >::type kron( V const& v, W const& w ) {
    return vector_kron_expression< V, W >( v, w ) ;
  }

} // namespace glas2

#endif
