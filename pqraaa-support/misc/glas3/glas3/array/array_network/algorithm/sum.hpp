//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_array_network_algorithm_sum_hpp
#define glas3_array_array_network_algorithm_sum_hpp

#include <glas3/concept/is.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/array_network/concept/tensor_network.hpp>

#include <glas3/array/array_network/container/tensor_network.hpp>

#include <type_traits>
#include <initializer_list>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

namespace glas3 {

template < typename T >
typename std::enable_if< is< TensorNetwork, T >::value, typename T::value_type >::type
sum( T const& t ) {
	auto inner_edges = t.inner_edges() ;
    for ( auto edge: t.outer_edges() ) {
    	inner_edges.push_back( edge ) ;
    }

    return tensor_network<typename T::array_type>( t.nodes2arrays(), t.nodes2edges(), inner_edges, {} )[0] ;
}

} // namespace glas3

#endif
