//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_backend_default_backend_array_dense_array_inner_prod_hpp
#define glas3_backend_default_backend_array_dense_array_inner_prod_hpp

#include <glas3/backend/default_backend/default_backend.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>
#include <glas3/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas3 {

template < typename X, typename Y >
typename std::enable_if< is< DenseArray, X >::value && is< DenseArray, Y >::value, decltype( typename X::value_type() * typename Y::value_type() ) >::type
inner_prod( default_backend, X const& x, Y const& y ) {
	assert( x.size() == y.size() ) ;
	decltype( typename X::value_type() * typename Y::value_type() ) sum = 0 ;
	for (typename X::size_type i = 0; i < x.size(); ++i) {
		sum += x[i] * y[i] ;
	}
	return sum ;
}

} // namespace glas3

#endif
