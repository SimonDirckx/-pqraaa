//  (C) Copyright Sam Corveleyn 2014.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_algorithm_full_hpp
#define glas3_array_dense_array_algorithm_full_hpp

#include <glas3/array/dense_array/concept/dense_array.hpp>
#include <glas3/concept/is.hpp>

#include <glas3/array/dense_array/container/dense_scalar.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_matrix.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <type_traits>

namespace glas3 {

template <typename X>
typename std::enable_if< is<DenseScalar, X>::value, dense_scalar<typename X::value_type> >::type
full( X const& x ) {
	return dense_scalar<typename X::value_type>( x ) ;
}

template <typename X>
typename std::enable_if< is<DenseVector, X>::value, dense_vector<typename X::value_type> >::type
full( X const& x ) {
	return dense_vector<typename X::value_type>( x ) ;
}

template <typename X>
typename std::enable_if< is<DenseMatrix, X>::value, dense_matrix<typename X::value_type> >::type
full( X const& x ) {
	return dense_matrix<typename X::value_type>( x ) ;
}

template <typename X>
typename std::enable_if< is<DenseArray, X>::value && ~is<DenseScalar, X>::value && ~is<DenseVector, X>::value && ~is<DenseMatrix, X>::value,
dense_array<typename X::value_type> >::type
full( X const& x ) {
	return dense_array<typename X::value_type>( x ) ;
}

} // namespace glas3

#endif
