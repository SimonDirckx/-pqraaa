//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_algorithm_squeeze_hpp
#define glas3_array_dense_array_algorithm_squeeze_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/type/reshape_index.hpp>
#include <glas3/array/type/indexed_array.hpp>

#include <type_traits>

namespace glas3 {

template < typename T >
typename std::enable_if< is< DenseArray, T >::value, indexed_array< T, reshape_index< typename T::size_type > > >::type
squeeze( T const& t ) {
	typename T::ndims_type ndims = 0 ;
	for ( typename T::ndims_type i = 0; i < t.shape().size(); ++i ) {
		if ( t.shape()[i] > 1 ) ndims += 1 ;
	}

	typename T::ndims_type idim = 0 ;
	dense_vector<typename T::size_type> shape( no_init(), ndims ) ;
	for ( typename T::ndims_type i = 0; i < t.shape().size(); ++i ) {
		if ( t.shape()[i] > 1 ) {
			shape[idim] = t.shape()[i] ;
			idim += 1 ;
		}
	}
	return indexed_array< T, reshape_index< typename T::size_type > >( t, reshape_index< typename T::size_type >( shape, t.shape() ) ) ;
}

} // namespace glas3

#endif
