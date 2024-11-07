//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_algorithm_tile_hpp
#define glas3_array_dense_array_algorithm_tile_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas3/array/dense_array/type/tile_index.hpp>
#include <glas3/array/type/indexed_array.hpp>

#include <initializer_list>
#include <type_traits>

namespace glas3 {

template < typename T, typename S >
typename std::enable_if< is< DenseArray, T >::value && is< Array, S >::value, indexed_array< T, tile_index< typename T::size_type > > >::type
tile( T const& t, S const& tile_shape ) { // remark: permute of empty array not possible
	return indexed_array< T, tile_index< typename T::size_type > >( t, tile_index< typename T::size_type >( t.shape(), tile_shape ) ) ;
}

template < typename T >
typename std::enable_if< is< DenseArray, T >::value, indexed_array< T, tile_index< typename T::size_type > > >::type
tile( T const& t, std::initializer_list< typename T::size_type > const& tile_shape ) { // remark: permute of empty array not possible
	return indexed_array< T, tile_index< typename T::size_type > >( t, tile_index< typename T::size_type >( t.shape(), tile_shape ) ) ;
}

} // namespace glas3

#endif
