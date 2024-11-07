//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_algorithm_vect_hpp
#define glas3_array_dense_array_algorithm_vect_hpp

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
vect( T const& t ) {
	return indexed_array< T, reshape_index< typename T::size_type > >( t, reshape_index< typename T::size_type >( t.size(), t.shape() ) ) ;
}

} // namespace glas3

#endif
