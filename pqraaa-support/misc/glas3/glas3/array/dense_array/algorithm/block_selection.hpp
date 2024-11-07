//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_algorithm_block_selection_hpp
#define glas3_array_dense_array_algorithm_block_selection_hpp

#include <glas3/concept/is.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas3/array/type/primitive_vector_wrapper.hpp>
#include <glas3/array/type/block_selection_array.hpp>

#include <initializer_list>
#include <type_traits>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>

namespace glas3 {

template < typename T >
typename std::enable_if< is< DenseArray, T >::value, block_selection_array< T > >::type
block_selection( T const& t, std::initializer_list< primitive_vector_wrapper< typename T::size_type > > const& s ) {
	assert( t.shape().size() == s.size() ) ;
    std::vector<boost::shared_ptr<primitive_vector_wrapper<typename T::size_type>>> selection_arrays( s.size() ) ;
	typename T::ndims_type k = 0 ;
	for ( auto s_it = s.begin(); s_it != s.end(); ++s_it, ++k ){
		selection_arrays[k] = boost::make_shared<primitive_vector_wrapper<typename T::size_type>>( s_it->shallow_copy() ) ;
	}
    return block_selection_array< T >( t, selection_arrays ) ;
}

template < typename T >
typename std::enable_if< is< DenseArray, T >::value, block_selection_array< T > >::type
block_selection( T const& t, std::vector< primitive_vector_wrapper< typename T::size_type > > const& s ) {
	assert( t.shape().size() == s.size() ) ;
    std::vector<boost::shared_ptr<primitive_vector_wrapper<typename T::size_type>>> selection_arrays( s.size() ) ;
	typename T::ndims_type k = 0 ;
	for ( auto s_it = s.begin(); s_it != s.end(); ++s_it, ++k ){
		selection_arrays[k] = boost::make_shared<primitive_vector_wrapper<typename T::size_type>>( s_it->shallow_copy() ) ;
	}
    return block_selection_array< T >( t, selection_arrays ) ;
}

} // namespace glas3

#endif
