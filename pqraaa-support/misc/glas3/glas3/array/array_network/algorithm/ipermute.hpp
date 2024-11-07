//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_array_network_algorithm_ipermute_hpp
#define glas3_array_array_network_algorithm_ipermute_hpp

#include <glas3/concept/is.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/array_network/concept/tensor_network.hpp>

#include <glas3/array/dense_array/container/dense_vector.hpp>

#include <glas3/array/array_network/algorithm/permute.hpp>

#include <type_traits>
#include <initializer_list>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

namespace glas3 {

template < typename T, typename D >
typename std::enable_if< is< TensorNetwork, T >::value && is< Array, D >::value, T >::type
ipermute( T const& t, D const& idim_order ) {
	dense_vector<typename T::ndims_type> dim_order( no_init(), idim_order.size() ) ;
	for ( std::ptrdiff_t i = 0; i < dim_order.size(); ++i ) {
		dim_order[idim_order[i]] = i ;
	}

	return permute( t, dim_order ) ;
}

template < typename T >
auto ipermute( T const& t, std::initializer_list< typename T::ndims_type > idim_order )
-> typename std::enable_if< is< TensorNetwork, T >::value, decltype( ipermute( t, dense_vector<typename T::ndims_type>( idim_order ) ) ) >::type {
	return ipermute( t, dense_vector<typename T::ndims_type>( idim_order ) ) ;
}

} // namespace glas3

#endif
