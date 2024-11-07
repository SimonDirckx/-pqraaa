//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_backend_default_backend_matrix_algorithm_norm_1_hpp
#define glas2_backend_default_backend_matrix_algorithm_norm_1_hpp

#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/backend/default_backend/vector/norm_1.hpp>
#include <glas2/backend/default_backend/default_backend.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cmath>

namespace glas2 {

  template <typename X>
  typename std::enable_if< is<DenseMatrix,X>::value
                         , decltype( std::abs(typename X::value_type()) )
                         >::type norm_1( default_backend, X const& x ) {
    decltype( std::abs(typename X::value_type()) ) mx = 0 ;
    for (typename X::size_type i=0; i<x.num_columns(); ++i) {
      mx = std::max( mx, glas2::norm_1( default_backend(), x( glas2::all(), i ) ) ) ;
    }
    return mx ;
  }
} // namespace glas2

#endif
