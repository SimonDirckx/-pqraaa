//  (C) Copyright Sam Corveleyn 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_algorithm_norm_hpp
#define glas3_array_dense_array_algorithm_norm_hpp

#include <glas3/array/algorithm/ndims.hpp>
#include <glas3/backend/current_backend.hpp>
#include <glas3/backend/default_backend/array/dense_array/norm.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>
#include <glas3/concept/is.hpp>
#include <type_traits>
#include <cmath>

namespace glas3 {

template <typename X, typename T>
auto norm( X const& x, T const& p ) -> typename std::enable_if< is<DenseArray, X>::value, decltype( norm( current_backend(), x, p ) ) >::type {
	return norm( current_backend(), x, p ) ;
}

template <typename X>
auto norm_0( X const& x ) -> typename std::enable_if< is<DenseArray, X>::value, decltype( norm_0( current_backend(), x ) ) >::type{
	return norm_0( current_backend(), x ) ;
}

template <typename X>
auto norm_1( X const& x ) -> typename std::enable_if< is<DenseArray, X>::value, decltype( norm_1( current_backend(), x ) ) >::type{
	return norm_1( current_backend(), x ) ;
}

template <typename X>
auto norm_2( X const& x ) -> typename std::enable_if< is<DenseArray, X>::value, decltype( norm_2( current_backend(), x ) ) >::type{
	return norm_2( current_backend(), x ) ;
}

template <typename X>
auto norm_inf( X const& x ) -> typename std::enable_if< is<DenseArray, X>::value, decltype( norm_inf( current_backend(), x ) ) >::type{
	return norm_inf( current_backend(), x ) ;
}

} // namespace glas3

#endif
