//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_array_network_algorithm_randomize_hpp
#define glas3_array_array_network_algorithm_randomize_hpp

#include <glas3/array/array_network/concept/tensor_network.hpp>
#include <glas3/concept/is.hpp>

#include <glas3/array/dense_array/algorithm/randomize.hpp>

#include <type_traits>

namespace glas3 {

  template <typename V, typename G>
  typename std::enable_if< is< TensorNetwork, V >::value >::type randomize( V const& v, G& g ) {
    for ( auto n: v.nodes2arrays() ) {
    	randomize( *n.second, g ) ;
    }
  }

} // namespace glas3

#endif
