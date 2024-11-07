//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_backend_default_backend_matrix_algorithm_norm_fro_squared_hpp
#define glas2_backend_default_backend_matrix_algorithm_norm_fro_squared_hpp

#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/backend/default_backend/default_backend.hpp>
#include <glas2/concept/is.hpp>
#include <glas2/concept/abs_squared.hpp>
#include <type_traits>
#include <cmath>

namespace glas2 {

  template <typename X>
  typename std::enable_if< is<DenseMatrix,X>::value //&& std::is_floating_point<typename X::value_type>::value
                         , decltype( std::abs(typename X::value_type()) )
                         >::type norm_fro_squared( default_backend, X const& x ) {
    decltype( std::abs(typename X::value_type()) ) sum = 0 ;
    for (typename X::size_type i=0; i<x.num_columns(); ++i) {
      for (typename X::size_type j=0; j<x.num_rows(); ++j) {
        sum += glas2::abs_squared( x(j,i) ) ;
      }
    }
    return sum ;
  }
} // namespace glas2

#endif
