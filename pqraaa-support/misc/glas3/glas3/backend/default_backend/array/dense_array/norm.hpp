//  (C) Copyright Sam Corveleyn 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_backend_default_backend_array_dense_array_algorithm_norm_hpp
#define glas3_backend_default_backend_array_dense_array_algorithm_norm_hpp

#include <glas3/backend/default_backend/default_backend.hpp>

#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas3/concept/is.hpp>

#include <type_traits>
#include <cmath>
#include <cassert>
#include <boost/type_traits/is_complex.hpp>
#include <complex>

namespace glas3 {

template <typename X>
typename std::enable_if< is<DenseArray, X>::value && std::is_arithmetic<typename X::value_type>::value, decltype( std::abs( typename X::value_type() ) ) >::type norm_0( default_backend, X const& x ) {
	typedef decltype( std::abs(typename X::value_type()) ) abs_type ;

	abs_type sum = 0 ;

	for (typename X::size_type i = 0; i < x.size(); ++i) { if ( x[i] != 0 ) sum += 1 ; }

	return sum ;
}

template <typename X>
typename std::enable_if< is<DenseArray, X>::value && boost::is_complex<typename X::value_type>::value, decltype( std::abs( typename X::value_type() ) ) >::type norm_0( default_backend, X const& x ) {
	typedef decltype( std::abs(typename X::value_type()) ) abs_type ;

	abs_type sum = 0 ;
	typename X::value_type complex_0(0, 0) ;

	for (typename X::size_type i = 0; i < x.size(); ++i) { if ( x[i] != complex_0 ) sum += 1 ; }

	return sum ;
}

template <typename X>
typename std::enable_if< is<DenseArray, X>::value, decltype( std::abs(typename X::value_type()) ) >::type norm_1( default_backend, X const& x ) {
	typedef decltype( std::abs(typename X::value_type()) ) abs_type ;

	abs_type sum = 0 ;

	for (typename X::size_type i = 0; i < x.size(); ++i) { sum += std::abs( x[i] ) ; }

	return sum ;
}

template <typename X>
typename std::enable_if< is<DenseArray, X>::value && std::is_arithmetic<typename X::value_type>::value, decltype( std::abs(typename X::value_type()) ) >::type norm_2( default_backend, X const& x ) {
	typedef decltype( std::abs(typename X::value_type()) ) abs_type ;

	abs_type sum = 0 ;

	for (typename X::size_type i = 0; i < x.size(); ++i) { sum += std::pow( x[i], 2 ) ; }

	return std::sqrt( sum ) ;
}

template <typename X>
typename std::enable_if< is<DenseArray, X>::value && boost::is_complex<typename X::value_type>::value, decltype( std::norm(typename X::value_type()) ) >::type norm_2( default_backend, X const& x ) {
	typedef decltype( std::norm(typename X::value_type()) ) abs_type ;

	abs_type sum = 0 ;

	for (typename X::size_type i = 0; i < x.size(); ++i) { sum += std::norm( x[i] ) ; }

	return std::sqrt( sum ) ;
}

template <typename X>
typename std::enable_if< is<DenseArray, X>::value, decltype( std::abs(typename X::value_type()) ) >::type norm_inf( default_backend, X const& x ) {
	typedef decltype( std::abs(typename X::value_type()) ) abs_type ;

	abs_type sum = 0 ;

	for (typename X::size_type i = 0; i < x.size(); ++i) { abs_type v( std::abs( x[i] ) ) ; if (v > sum) sum = v ; }

	return sum ;
}

template <typename X, typename T>
typename std::enable_if< is<DenseArray, X>::value, decltype( std::abs(typename X::value_type()) ) >::type norm( default_backend, X const& x, T const& p ) {
    typedef decltype( std::abs(typename X::value_type()) ) abs_type ;
    assert( p >= 0 ) ;

    abs_type sum ;

    if ( p == 0 ) { sum = norm_0( x ) ; }
    else if ( p == 1 ) { sum = norm_1( x ) ; }
    else if ( p == 2 ) { sum = norm_2( x ) ; }
    else if ( std::isinf(p) ) { sum = norm_inf( x ) ; }
    else {
    	sum = 0 ;
        for (typename X::size_type i = 0; i < x.size(); ++i) { sum += std::pow( std::abs( x[i] ), p ) ; } // could be more efficient for complex numbers
        sum = std::pow ( sum,  1. / p ) ;
    }

    return sum ;
}

} // namespace glas3

#endif
