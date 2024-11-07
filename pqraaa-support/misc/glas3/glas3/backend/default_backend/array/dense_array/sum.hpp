//  (C) Copyright Sam Corveleyn 2014.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_backend_default_backend_array_dense_array_sum_hpp
#define glas3_backend_default_backend_array_dense_array_sum_hpp

#include <glas3/backend/default_backend/default_backend.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>
#include <glas3/concept/is.hpp>
#include <type_traits>
#include <cmath>

namespace glas3 {

template <typename X>
typename std::enable_if< is<DenseArray, X>::value, typename X::value_type>::type sum( default_backend, X const& x ) {
	typename X::value_type sum = 0 ;

	for (typename X::size_type i = 0; i < x.size(); ++i) { sum += x[i] ; }

	return sum ;
}

} // namespace glas3

#endif
