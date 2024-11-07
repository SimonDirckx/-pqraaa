//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_backend_blas_backend_vector_algorithm_sum_hpp
#define glas2_backend_blas_backend_vector_algorithm_sum_hpp

#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/backend/blas_backend/blas_backend.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cmath>

namespace glas2 {

  template <typename X>
  typename std::enable_if< is<DenseVector,X>::value
                         , typename X::value_type
                         >::type sum( blas_backend, X const& x ) {
    typename X::value_type sum = 0 ;
    for (typename X::size_type i=0; i<x.size(); ++i) {
      sum += x(i) ;
    }
    return sum ;
  }

} // namespace glas2

#endif
