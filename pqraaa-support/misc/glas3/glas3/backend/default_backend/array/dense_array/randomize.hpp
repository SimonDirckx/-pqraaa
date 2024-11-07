//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_backend_default_backend_array_dense_array_randomize_hpp
#define glas3_backend_default_backend_array_dense_array_randomize_hpp

#include <glas3/backend/default_backend/default_backend.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>
#include <glas3/concept/is.hpp>

#include <type_traits>
#include <boost/type_traits/is_complex.hpp>

namespace glas3 {

template <typename V, typename G>
typename std::enable_if< is< DenseArray, V >::value && std::is_arithmetic< typename V::value_type >::value >::type randomize( default_backend, V const& v, G& g ) {
	for ( typename V::size_type i = 0; i < v.size(); ++i) {
		v[i] = g() ;
	}
}

template <typename V, typename G>
typename std::enable_if< is< DenseArray, V >::value && boost::is_complex< typename V::value_type >::value >::type randomize( default_backend, V const& v, G& g ) {
	for ( typename V::size_type i = 0; i < v.size(); ++i) {
		v[i].real( g() ) ;
		v[i].imag( g() ) ;
	}
}

} // namespace glas3

#endif
