//  (C) Copyright Sam Corveleyn 2014.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_algorithm_find_hpp
#define glas3_array_dense_array_algorithm_find_hpp

#include <glas3/concept/is.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas3/array/dense_array/container/dense_vector.hpp>

#include <vector>
#include <type_traits>

namespace glas3 {

template <typename X>
typename std::enable_if< is< DenseArray, X >::value, dense_vector<typename X::size_type> >::type
find( X const& x ) { // needs more efficient implementation
	std::vector<typename X::size_type> ind ;
	for ( typename X::size_type i = 0; i < x.size(); ++i ) {
		if ( x[i] != 0 ) ind.push_back( i ) ;
	}

	dense_vector<typename X::size_type> v( glas3::no_init(), ind.size() ) ;
	for ( std::ptrdiff_t i = 0; i < ind.size(); ++i ) {
		v[i] = ind[i] ;
	}

	return v ;
}

} // namespace glas3

#endif



