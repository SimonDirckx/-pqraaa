//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_algorithm_transpose_hpp
#define glas3_array_dense_array_algorithm_transpose_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas3/array/dense_array/container/dense_vector.hpp>

#include <glas3/array/dense_array/type/range.hpp>

#include <glas3/array/dense_array/algorithm/permute.hpp>

#include <type_traits>

namespace glas3 {

template < typename T >
auto transpose( T const& t ) -> typename std::enable_if< is< DenseArray, T >::value, decltype( permute( t, glas3::range<typename T::ndims_type>() ) ) >::type {
	return permute( t, glas3::range<typename T::ndims_type>( t.shape().size() - 1, -1, -1 ) ) ;
}

} // namespace glas3

#endif
