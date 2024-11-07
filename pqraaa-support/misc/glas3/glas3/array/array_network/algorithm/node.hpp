//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_array_network_algorithm_node_hpp
#define glas3_array_array_network_algorithm_node_hpp

#include <glas3/concept/is.hpp>

#include <glas3/array/array_network/concept/tensor_network.hpp>

#include <type_traits>

namespace glas3 {

  template <typename X>
  auto node( X const& x, typename X::ndims_type const& k ) -> typename std::enable_if< is<TensorNetwork, X>::value, decltype( x.nodes2arrays().at( k )->shallow_copy() ) >::type {
    return x.nodes2arrays().at( k )->shallow_copy() ; ;
  }

} // namespace glas3

#endif
