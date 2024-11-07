//  (C) Copyright Sam Corveleyn 2014.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_algorithm_sum_hpp
#define glas3_array_dense_array_algorithm_sum_hpp

#include <glas3/backend/current_backend.hpp>
#include <glas3/backend/default_backend/array/dense_array/sum.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>
#include <glas3/concept/is.hpp>
#include <type_traits>
#include <cmath>

namespace glas3 {

template <typename X>
auto sum( X const& x ) -> typename std::enable_if< is<DenseArray, X>::value, decltype( sum( current_backend(), x ) ) >::type{
	return sum( current_backend(), x ) ;
}

} // namespace glas3

#endif



