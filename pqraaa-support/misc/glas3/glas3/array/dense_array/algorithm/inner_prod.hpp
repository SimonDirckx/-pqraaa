//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_algorithm_inner_prod_hpp
#define glas3_array_dense_array_algorithm_inner_prod_hpp

#include <glas3/backend/current_backend.hpp>
#include <glas3/backend/default_backend/array/dense_array/inner_prod.hpp>

#include <glas3/array/dense_array/concept/dense_array.hpp>
#include <glas3/concept/is.hpp>

#include <type_traits>
#include <cmath>

namespace glas3 {

template <typename X, typename Y>
auto inner_prod( X const& x, Y const& y ) -> typename std::enable_if< is<DenseArray, X>::value && is<DenseArray, Y>::value, decltype( inner_prod( current_backend(), x, y ) ) >::type {
	return inner_prod( current_backend(), x, y ) ;
}

} // namespace glas3

#endif
