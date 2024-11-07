//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_array_network_algorithm_rank_hpp
#define glas3_array_array_network_algorithm_rank_hpp

#include <glas3/concept/is.hpp>

#include <glas3/array/array_network/concept/tensor_network.hpp>

#include <type_traits>

namespace glas3 {

  template <typename X>
  auto rank( X const& x ) -> typename std::enable_if< is<TensorNetwork, X>::value, decltype( x.rank().shallow_copy() ) >::type {
    return x.rank().shallow_copy() ;
  }

} // namespace glas3

#endif
