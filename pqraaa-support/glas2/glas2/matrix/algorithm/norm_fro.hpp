//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_algorithm_norm_fro_hpp
#define glas2_matrix_algorithm_norm_fro_hpp

#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/matrix/algorithm/norm_fro_squared.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cmath>

namespace glas2 {

  template <typename X>
  typename std::enable_if< is<DenseMatrix,X>::value
                         , decltype( std::abs(typename X::value_type()) )
                         >::type norm_fro( X const& x ) {
    return std::sqrt( norm_fro_squared( x ) ) ;
  }

  template <typename B, typename X>
  typename std::enable_if< is<DenseMatrix,X>::value
                         , decltype( std::abs(typename X::value_type()) )
                         >::type norm_fro( B const& backend, X const& x ) {
    return std::sqrt( norm_fro_squared( backend, x ) ) ;
  }

} // namespace glas2

#endif
