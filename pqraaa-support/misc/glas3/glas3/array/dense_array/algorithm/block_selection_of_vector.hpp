//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_algorithm_block_selection_of_vector_hpp
#define glas3_array_dense_array_algorithm_block_selection_of_vector_hpp

#include <glas3/concept/is.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas3/array/type/primitive_vector_wrapper.hpp>

#include <glas3/array/type/block_selection_vector.hpp>

#include <type_traits>

namespace glas3 {

template < typename T >
typename std::enable_if< is< DenseVector, T >::value, block_selection_vector< T > >::type
block_selection_of_vector( T const& t, primitive_vector_wrapper<typename T::size_type> const& s ) {
	return block_selection_vector< T >( t, s ) ;
}

} // namespace glas3

#endif
