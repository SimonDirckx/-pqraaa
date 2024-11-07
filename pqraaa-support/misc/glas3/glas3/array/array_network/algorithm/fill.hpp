//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_array_network_algorithm_fill_hpp
#define glas3_array_array_network_algorithm_fill_hpp

#include <glas3/concept/is.hpp>
#include <glas3/array/array_network/concept/tensor_network.hpp>

#include <glas3/array/dense_array/algorithm/fill.hpp>

#include <type_traits>

namespace glas3 {

template <typename T, typename V>
typename std::enable_if< is< TensorNetwork, T >::value >::type fill( T const& t, V const& v ) {
    for ( auto n: t.nodes2arrays() ) {
    	fill( *n.second, v ) ;
    }
}

} // namespace glas3

#endif
