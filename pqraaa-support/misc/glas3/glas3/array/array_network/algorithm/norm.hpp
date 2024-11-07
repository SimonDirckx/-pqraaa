//  (C) Copyright Sam Corveleyn 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_array_network_algorithm_norm_hpp
#define glas3_array_array_network_algorithm_norm_hpp

#include <glas3/concept/is.hpp>
#include <glas3/array/array_network/concept/tensor_network.hpp>

#include <glas3/array/array_network/algorithm/inner_prod.hpp>

#include <cmath>
#include <type_traits>

namespace glas3 {

  template <typename X, typename T>
  auto norm( X const& x, T const& p ) -> typename std::enable_if< is<TensorNetwork, X>::value, decltype( std::sqrt( inner_prod( x, x ) ) ) >::type {
	    assert( p == 2 ) ;

	    return std::sqrt( inner_prod( x, x ) ) ;
  }

  template <typename X>
  auto norm_2( X const& x ) -> typename std::enable_if< is<TensorNetwork, X>::value, decltype( std::sqrt( inner_prod( x, x ) ) ) >::type{
	  return std::sqrt( inner_prod( x, x ) ) ;
  }

} // namespace glas3

#endif
