//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_algorithm_reshape_hpp
#define glas3_array_dense_array_algorithm_reshape_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/type/reshape_index.hpp>
#include <glas3/array/type/indexed_array.hpp>

#include <initializer_list>
#include <type_traits>

namespace glas3 {

template < typename T, typename S >
typename std::enable_if< is< DenseArray, T >::value && is< Array, S >::value, indexed_array< T, reshape_index< typename T::size_type > > >::type
reshape( T const& t, S const& s ) {

	typename S::value_type size_new = 1;
	for ( typename S::size_type k = 0; k < s.size(); ++k ) { size_new *= s[k] ; }
	assert( size_new == t.size() ) ; // remark: reshape of empty array not possible

	return indexed_array< T, reshape_index< typename T::size_type > >( t, reshape_index< typename T::size_type >( s, t.shape() ) ) ;
}

template < typename T >
typename std::enable_if< is< DenseArray, T >::value, indexed_array< T, reshape_index< typename T::size_type > > >::type
reshape( T const& t, std::initializer_list< typename T::size_type > const& s ) {

	typename T::size_type size_new = 1;
	auto s_ptr = s.begin() ;
	for ( std::ptrdiff_t k = 0; k < s.size(); ++k, ++s_ptr ) { size_new *= *s_ptr ; }
	assert( size_new == t.size() ) ; // remark: reshape of empty array not possible

	return indexed_array< T, reshape_index< typename T::size_type > >( t, reshape_index< typename T::size_type >( s, t.shape() ) ) ;
}

} // namespace glas3

#endif
