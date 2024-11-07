//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_algorithm_norm_2_hpp
#define glas2_vector_algorithm_norm_2_hpp

#include <glas2/vector/algorithm/norm_2_squared.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cmath>

namespace glas2 {

  template <typename X>
  typename std::enable_if< is<DenseVector,X>::value
                         , decltype( std::abs(typename X::value_type()) )
                         >::type norm_2( X const& x ) {
    return norm_2( current_backend(), x ) ;
  }

  template <typename B, typename X>
  typename std::enable_if< is<DenseVector,X>::value
                         , decltype( std::abs(typename X::value_type()) )
                         >::type norm_2( B const& backend, X const& x ) {
    return std::sqrt( norm_2_squared( backend, x ) ) ;
  }

} // namespace glas2

#endif
