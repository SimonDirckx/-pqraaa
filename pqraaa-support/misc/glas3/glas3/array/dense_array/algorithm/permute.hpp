//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_algorithm_permute_hpp
#define glas3_array_dense_array_algorithm_permute_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas3/array/dense_array/type/permute_index.hpp>
#include <glas3/array/type/indexed_array.hpp>

#include <type_traits>
#include <initializer_list>

namespace glas3 {

template < typename T, typename D >
typename std::enable_if< is< DenseArray, T >::value && is< Array, D >::value, indexed_array< T, permute_index< typename T::size_type > > >::type
permute( T const& t, D const& dim_order ) { // remark: permute of empty array not possible
    assert( t.shape().size() == dim_order.size() ) ;
	return indexed_array< T, permute_index< typename T::size_type > >( t, permute_index< typename T::size_type >( t.shape(), dim_order ) ) ;
}

template < typename T >
typename std::enable_if< is< DenseArray, T >::value, indexed_array< T, permute_index< typename T::size_type > > >::type
permute( T const& t, std::initializer_list< typename T::size_type > const& dim_order ) { // remark: permute of empty array not possible
	assert( t.shape().size() == dim_order.size() ) ;
	return indexed_array< T, permute_index< typename T::size_type > >( t, permute_index< typename T::size_type >( t.shape(), dim_order ) ) ;
}

} // namespace glas3

#endif
