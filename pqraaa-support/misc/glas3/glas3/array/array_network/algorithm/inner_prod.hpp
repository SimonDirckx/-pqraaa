//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_array_network_algorithm_inner_prod_hpp
#define glas3_array_array_network_algorithm_inner_prod_hpp

#include <glas3/concept/is.hpp>
#include <glas3/array/array_network/concept/tensor_network.hpp>

#include <glas3/array/array_network/algorithm/ttt.hpp>

#include <glas3/array/dense_array/type/range.hpp>

#include <type_traits>

namespace glas3 {

  template <typename X, typename Y>
  auto inner_prod( X const& x, Y const& y )
  -> typename std::enable_if< is<TensorNetwork, X>::value && is<TensorNetwork, Y>::value, decltype( ttt( x, y, range<typename X::ndims_type>( 0, x.shape().size() ), range<typename Y::ndims_type>( 0, y.shape().size() ) )[0] ) >::type {
      return ttt( x, y, range<typename X::ndims_type>( 0, x.shape().size() ), range<typename Y::ndims_type>( 0, y.shape().size() ) )[0] ;
  }

} // namespace glas3

#endif
