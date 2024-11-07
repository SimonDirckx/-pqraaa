//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_algorithm_ipermute_hpp
#define glas3_array_dense_array_algorithm_ipermute_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/type/vector_wrapper.hpp>
#include <glas3/array/dense_array/type/permute_index.hpp>
#include <glas3/array/type/indexed_array.hpp>

#include <type_traits>

namespace glas3 {

template < typename T >
typename std::enable_if< is< DenseArray, T >::value, indexed_array< T, permute_index< typename T::size_type > > >::type
ipermute( T const& t, vector_wrapper<typename T::ndims_type> const& idim_order ) { // remark: permute of empty array not possible
	dense_vector<typename T::ndims_type> dim_order( no_init(), idim_order.size() ) ;
	for ( std::ptrdiff_t i = 0; i < dim_order.size(); ++i ) {
		dim_order[idim_order[i]] = i ;
	}
	return indexed_array< T, permute_index< typename T::size_type > >( t, permute_index< typename T::size_type >( t.shape(), dim_order ) ) ;
}

} // namespace glas3

#endif
