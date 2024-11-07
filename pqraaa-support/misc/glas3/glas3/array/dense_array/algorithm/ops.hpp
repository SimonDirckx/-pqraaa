//  (C) Copyright Sam Corveleyn 2014.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_algorithm_ops_hpp
#define glas3_array_dense_array_algorithm_ops_hpp

#include <glas3/array/dense_array/type/entrywise_operation.hpp>
#include <glas3/array/dense_array/type/constant_array.hpp>

#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas3/concept/is.hpp>
#include <type_traits>
#include <cmath>
#include <functional>
#include <boost/type_traits/is_complex.hpp>
#include <complex>

namespace glas3 {

template <typename Op, typename A1>
typename std::enable_if< is<DenseArray, A1>::value, unary_operation< Op, A1 > >::type
for_each( Op const& op, A1 const& a1 ) {
    return unary_operation< Op, A1 >( op, a1 ) ;
}

template <typename Op, typename A1, typename A2>
typename std::enable_if< is<DenseArray, A1>::value && is<DenseArray, A2>::value, binary_operation< Op, A1, A2 > >::type
for_each( Op const& op, A1 const& a1, A2 const& a2 ) {
    return binary_operation< Op, A1, A2 >( op, a1, a2 ) ;
}

template <typename X, typename Y>
struct ops_multiply {
	typedef decltype( X() * Y() ) result_type ;
	result_type operator() ( X const& x, Y const& y ) const { return x * y ; }
} ;

template <typename A1, typename A2>
typename std::enable_if< is<DenseArray, A1>::value && is<DenseArray, A2>::value,
binary_operation< ops_multiply<typename A1::value_type, typename A2::value_type>, A1, A2 > >::type
operator*( A1 const& a1, A2 const& a2 ) {
	return binary_operation< ops_multiply<typename A1::value_type, typename A2::value_type>, A1, A2 > ( ops_multiply<typename A1::value_type, typename A2::value_type>(), a1, a2 ) ;
}

template <typename X, typename Y>
struct ops_divide {
	typedef decltype( X() / Y() ) result_type ;
	result_type operator() ( X const& x, Y const& y ) const { return x / y ; }
} ;

template <typename A1, typename A2>
typename std::enable_if< is<DenseArray, A1>::value && is<DenseArray, A2>::value,
binary_operation< ops_divide<typename A1::value_type, typename A2::value_type>, A1, A2 > >::type
operator/( A1 const& a1, A2 const& a2 ) {
	return binary_operation< ops_divide<typename A1::value_type, typename A2::value_type>, A1, A2 > ( ops_divide<typename A1::value_type, typename A2::value_type>(), a1, a2 ) ;
}

template <typename X, typename Y>
struct ops_plus {
	typedef decltype( X() + Y() ) result_type ;
	result_type operator() ( X const& x, Y const& y ) const { return x + y ; }
} ;

template <typename A1, typename A2>
typename std::enable_if< is<DenseArray, A1>::value && is<DenseArray, A2>::value,
binary_operation< ops_plus<typename A1::value_type, typename A2::value_type>, A1, A2 > >::type
operator+( A1 const& a1, A2 const& a2 ) {
	return binary_operation< ops_plus<typename A1::value_type, typename A2::value_type>, A1, A2 > ( ops_plus<typename A1::value_type, typename A2::value_type>(), a1, a2 ) ;
}

template <typename X, typename Y>
struct ops_minus {
	typedef decltype( X() - Y() ) result_type ;
	result_type operator() ( X const& x, Y const& y ) const { return x - y ; }
} ;

template <typename A1, typename A2>
typename std::enable_if< is<DenseArray, A1>::value && is<DenseArray, A2>::value,
binary_operation< ops_minus<typename A1::value_type, typename A2::value_type>, A1, A2 > >::type
operator-( A1 const& a1, A2 const& a2 ) {
	return binary_operation< ops_minus<typename A1::value_type, typename A2::value_type>, A1, A2 > ( ops_minus<typename A1::value_type, typename A2::value_type>(), a1, a2 ) ;
}

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value, A1 >::type
operator+( A1 const& a1 ) {
	return a1.shallow_copy() ;
}

template <typename Y>
struct ops_uminus {
	typedef decltype( - Y() ) result_type ;
	result_type operator() ( Y const& y ) const { return - y ; }
} ;

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value, unary_operation< ops_uminus<typename A1::value_type>, A1 > >::type
operator-( A1 const& a1 ) {
	return unary_operation< ops_uminus<typename A1::value_type>, A1 >( ops_uminus<typename A1::value_type>(), a1 )  ;
}

template <typename X, typename Y>
struct ops_mod {
	typedef decltype( X() % Y() ) result_type ;
	result_type operator() ( X const& x, Y const& y ) const { return x % y ; }
} ;

template <typename A1, typename A2>
typename std::enable_if< is<DenseArray, A1>::value && is<DenseArray, A2>::value,
binary_operation< ops_mod<typename A1::value_type, typename A2::value_type>, A1, A2 > >::type
operator%( A1 const& a1, A2 const& a2 ) {
	return binary_operation< ops_mod<typename A1::value_type, typename A2::value_type>, A1, A2 > ( ops_mod<typename A1::value_type, typename A2::value_type>(), a1, a2 ) ;
}

template <typename X, typename Y>
struct ops_equal {
	typedef decltype( X() == Y() ) result_type ;
	result_type operator() ( X const& x, Y const& y ) const { return x == y ; }
} ;

template <typename A1, typename A2>
typename std::enable_if< is<DenseArray, A1>::value && is<DenseArray, A2>::value,
binary_operation< ops_equal<typename A1::value_type, typename A2::value_type>, A1, A2 > >::type
operator==( A1 const& a1, A2 const& a2 ) {
	return binary_operation< ops_equal<typename A1::value_type, typename A2::value_type>, A1, A2 > ( ops_equal<typename A1::value_type, typename A2::value_type>(), a1, a2 ) ;
}

template <typename X, typename Y>
struct ops_not_equal {
	typedef decltype( X() != Y() ) result_type ;
	result_type operator() ( X const& x, Y const& y ) const { return x != y ; }
} ;

template <typename A1, typename A2>
typename std::enable_if< is<DenseArray, A1>::value && is<DenseArray, A2>::value,
binary_operation< ops_not_equal<typename A1::value_type, typename A2::value_type>, A1, A2 > >::type
operator!=( A1 const& a1, A2 const& a2 ) {
	return binary_operation< ops_not_equal<typename A1::value_type, typename A2::value_type>, A1, A2 > ( ops_not_equal<typename A1::value_type, typename A2::value_type>(), a1, a2 ) ;
}

template <typename X, typename Y>
struct ops_greater_than {
	typedef decltype( X() > Y() ) result_type ;
	result_type operator() ( X const& x, Y const& y ) const { return x > y ; }
} ;

template <typename A1, typename A2>
typename std::enable_if< is<DenseArray, A1>::value && is<DenseArray, A2>::value,
binary_operation< ops_greater_than<typename A1::value_type, typename A2::value_type>, A1, A2 > >::type
operator>( A1 const& a1, A2 const& a2 ) {
	return binary_operation< ops_greater_than<typename A1::value_type, typename A2::value_type>, A1, A2 > ( ops_greater_than<typename A1::value_type, typename A2::value_type>(), a1, a2 ) ;
}

template <typename X, typename Y>
struct ops_greater_than_or_equal {
	typedef decltype( X() >= Y() ) result_type ;
	result_type operator() ( X const& x, Y const& y ) const { return x >= y ; }
} ;

template <typename A1, typename A2>
typename std::enable_if< is<DenseArray, A1>::value && is<DenseArray, A2>::value,
binary_operation< ops_greater_than_or_equal<typename A1::value_type, typename A2::value_type>, A1, A2 > >::type
operator>=( A1 const& a1, A2 const& a2 ) {
	return binary_operation< ops_greater_than_or_equal<typename A1::value_type, typename A2::value_type>, A1, A2 > ( ops_greater_than_or_equal<typename A1::value_type, typename A2::value_type>(), a1, a2 ) ;
}

template <typename X, typename Y>
struct ops_smaller_than {
	typedef decltype( X() < Y() ) result_type ;
	result_type operator() ( X const& x, Y const& y ) const { return x < y ; }
} ;

template <typename A1, typename A2>
typename std::enable_if< is<DenseArray, A1>::value && is<DenseArray, A2>::value,
binary_operation< ops_smaller_than<typename A1::value_type, typename A2::value_type>, A1, A2 > >::type
operator<( A1 const& a1, A2 const& a2 ) {
	return binary_operation< ops_smaller_than<typename A1::value_type, typename A2::value_type>, A1, A2 > ( ops_smaller_than<typename A1::value_type, typename A2::value_type>(), a1, a2 ) ;
}

template <typename X, typename Y>
struct ops_smaller_than_or_equel {
	typedef decltype( X() <= Y() ) result_type ;
	result_type operator() ( X const& x, Y const& y ) const { return x <= y ; }
} ;

template <typename A1, typename A2>
typename std::enable_if< is<DenseArray, A1>::value && is<DenseArray, A2>::value,
binary_operation< ops_smaller_than_or_equel<typename A1::value_type, typename A2::value_type>, A1, A2 > >::type
operator<=( A1 const& a1, A2 const& a2 ) {
	return binary_operation< ops_smaller_than_or_equel<typename A1::value_type, typename A2::value_type>, A1, A2 > ( ops_smaller_than_or_equel<typename A1::value_type, typename A2::value_type>(), a1, a2 ) ;
}

template <typename Y>
struct ops_not {
	typedef decltype( ! Y() ) result_type ;
	result_type operator() ( Y const& y ) const { return ! y ; }
} ;

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value, unary_operation< ops_not<typename A1::value_type>, A1 > >::type
operator!( A1 const& a1 ) {
	return unary_operation< ops_not<typename A1::value_type>, A1 >( ops_not<typename A1::value_type>(), a1 )  ;
}

template <typename X, typename Y>
struct ops_and {
	typedef decltype( X() && Y() ) result_type ;
	result_type operator() ( X const& x, Y const& y ) const { return x && y ; }
} ;

template <typename A1, typename A2>
typename std::enable_if< is<DenseArray, A1>::value && is<DenseArray, A2>::value,
binary_operation< ops_and<typename A1::value_type, typename A2::value_type>, A1, A2 > >::type
operator&&( A1 const& a1, A2 const& a2 ) {
	return binary_operation< ops_and<typename A1::value_type, typename A2::value_type>, A1, A2 > ( ops_and<typename A1::value_type, typename A2::value_type>(), a1, a2 ) ;
}

template <typename X, typename Y>
struct ops_or {
	typedef decltype( X() || Y() ) result_type ;
	result_type operator() ( X const& x, Y const& y ) const { return x || y ; }
} ;

template <typename A1, typename A2>
typename std::enable_if< is<DenseArray, A1>::value && is<DenseArray, A2>::value,
binary_operation< ops_or<typename A1::value_type, typename A2::value_type>, A1, A2 > >::type
operator||( A1 const& a1, A2 const& a2 ) {
	return binary_operation< ops_or<typename A1::value_type, typename A2::value_type>, A1, A2 > ( ops_or<typename A1::value_type, typename A2::value_type>(), a1, a2 ) ;
}

template <typename X, typename Y>
struct ops_pow {
	typedef decltype( std::pow( X(), Y() ) ) result_type ;
	result_type operator() ( X const& x, Y const& y ) const { return std::pow( x, y ) ; }
} ;

template <typename A1, typename A2>
typename std::enable_if< is<DenseArray, A1>::value && is<DenseArray, A2>::value,
binary_operation< ops_pow<typename A1::value_type, typename A2::value_type>, A1, A2 > >::type
pow( A1 const& a1, A2 const& a2 ) {
	return binary_operation< ops_pow<typename A1::value_type, typename A2::value_type>, A1, A2 > ( ops_pow<typename A1::value_type, typename A2::value_type>(), a1, a2 ) ;
}

template <typename Y>
struct ops_exp {
	typedef decltype( std::exp( Y() ) ) result_type ;
	result_type operator() ( Y const& y ) const { return std::exp( y ) ; }
} ;

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value, unary_operation< ops_exp<typename A1::value_type>, A1 > >::type
exp( A1 const& a1 ) {
	return unary_operation< ops_exp<typename A1::value_type>, A1 >( ops_exp<typename A1::value_type>(), a1 )  ;
}

template <typename Y>
struct ops_log {
	typedef decltype( std::log( Y() ) ) result_type ;
	result_type operator() ( Y const& y ) const { return std::log( y ) ; }
} ;

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value, unary_operation< ops_log<typename A1::value_type>, A1 > >::type
log( A1 const& a1 ) {
	return unary_operation< ops_log<typename A1::value_type>, A1 >( ops_log<typename A1::value_type>(), a1 )  ;
}

template <typename Y>
struct ops_log10 {
	typedef decltype( std::log10( Y() ) ) result_type ;
	result_type operator() ( Y const& y ) const { return std::log10( y ) ; }
} ;

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value, unary_operation< ops_log10<typename A1::value_type>, A1 > >::type
log10( A1 const& a1 ) {
	return unary_operation< ops_log10<typename A1::value_type>, A1 >( ops_log10<typename A1::value_type>(), a1 )  ;
}

template <typename Y>
struct ops_sqrt {
	typedef decltype( std::sqrt( Y() ) ) result_type ;
	result_type operator() ( Y const& y ) const { return std::sqrt( y ) ; }
} ;

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value, unary_operation< ops_sqrt<typename A1::value_type>, A1 > >::type
sqrt( A1 const& a1 ) {
	return unary_operation< ops_sqrt<typename A1::value_type>, A1 >( ops_sqrt<typename A1::value_type>(), a1 )  ;
}

template <typename Y>
struct ops_conj {
	typedef decltype( std::conj( Y() ) ) result_type ;
	result_type operator() ( Y const& y ) const { return std::conj( y ) ; }
} ;

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value && std::is_arithmetic<typename A1::value_type>::value, A1 >::type
conj( A1 const& a1 ) {
	return a1.shallow_copy()  ;
}

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value && boost::is_complex<typename A1::value_type>::value, unary_operation< ops_conj<typename A1::value_type>, A1 > >::type
conj( A1 const& a1 ) {
	return unary_operation< ops_conj<typename A1::value_type>, A1 >( ops_conj<typename A1::value_type>(), a1 )  ;
}

template <typename Y>
struct ops_imag {
	typedef decltype( std::imag( Y() ) ) result_type ;
	result_type operator() ( Y const& y ) const { return std::imag( y ) ; }
} ;

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value
                    && ! is<DenseMatrix, A1>::value
                    && ! is<DenseVector, A1>::value
                    && ! is<DenseScalar, A1>::value && std::is_arithmetic<typename A1::value_type>::value, constant_array<typename A1::value_type> >::type
imag( A1 const& a1 ) {
	return constant_array<typename A1::value_type>( 0, a1.shape() ) ;
}

template <typename A1>
typename std::enable_if< is<DenseScalar, A1>::value && std::is_arithmetic<typename A1::value_type>::value, constant_array<typename A1::value_type, 0> >::type
imag( A1 const& a1 ) {
	return constant_array<typename A1::value_type, 0>( 0 ) ;
}

template <typename A1>
typename std::enable_if< is<DenseVector, A1>::value && std::is_arithmetic<typename A1::value_type>::value, constant_array<typename A1::value_type, 1> >::type
imag( A1 const& a1 ) {
	return constant_array<typename A1::value_type, 1>( 0, a1.shape() ) ;
}

template <typename A1>
typename std::enable_if< is<DenseMatrix, A1>::value && std::is_arithmetic<typename A1::value_type>::value, constant_array<typename A1::value_type, 2> >::type
imag( A1 const& a1 ) {
	return constant_array<typename A1::value_type, 2>( 0, a1.shape() ) ;
}

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value && boost::is_complex<typename A1::value_type>::value, unary_operation< ops_imag<typename A1::value_type>, A1 > >::type
imag( A1 const& a1 ) {
	return unary_operation< ops_imag<typename A1::value_type>, A1 >( ops_imag<typename A1::value_type>(), a1 )  ;
}

template <typename Y>
struct ops_real {
	typedef decltype( std::real( Y() ) ) result_type ;
	result_type operator() ( Y const& y ) const { return std::real( y ) ; }
} ;

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value && std::is_arithmetic<typename A1::value_type>::value, A1 >::type
real( A1 const& a1 ) {
	return a1.shallow_copy()  ;
}

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value && boost::is_complex<typename A1::value_type>::value, unary_operation< ops_real<typename A1::value_type>, A1 > >::type
real( A1 const& a1 ) {
	return unary_operation< ops_real<typename A1::value_type>, A1 >( ops_real<typename A1::value_type>(), a1 )  ;
}

template <typename Y>
struct ops_abs {
	typedef decltype( std::abs( Y() ) ) result_type ;
	result_type operator() ( Y const& y ) const { return std::abs( y ) ; }
} ;

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value, unary_operation< ops_abs<typename A1::value_type>, A1 > >::type
abs( A1 const& a1 ) {
	return unary_operation< ops_abs<typename A1::value_type>, A1 >( ops_abs<typename A1::value_type>(), a1 )  ;
}

template <typename Y>
struct ops_arg {
	typedef decltype( std::arg( Y() ) ) result_type ;
	result_type operator() ( Y const& y ) const { return std::arg( y ) ; }
} ;

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value
                    && ! is<DenseMatrix, A1>::value
                    && ! is<DenseVector, A1>::value
                    && ! is<DenseScalar, A1>::value && std::is_arithmetic<typename A1::value_type>::value, constant_array<typename A1::value_type> >::type
arg( A1 const& a1 ) {
	return constant_array<typename A1::value_type>( 0, a1.shape() ) ;
}

template <typename A1>
typename std::enable_if< is<DenseScalar, A1>::value && std::is_arithmetic<typename A1::value_type>::value, constant_array<typename A1::value_type, 0> >::type
arg( A1 const& a1 ) {
	return constant_array<typename A1::value_type, 0>( 0 ) ;
}

template <typename A1>
typename std::enable_if< is<DenseVector, A1>::value && std::is_arithmetic<typename A1::value_type>::value, constant_array<typename A1::value_type, 1> >::type
arg( A1 const& a1 ) {
	return constant_array<typename A1::value_type, 1>( 0, a1.shape() ) ;
}

template <typename A1>
typename std::enable_if< is<DenseMatrix, A1>::value && std::is_arithmetic<typename A1::value_type>::value, constant_array<typename A1::value_type, 2> >::type
arg( A1 const& a1 ) {
	return constant_array<typename A1::value_type, 2>( 0, a1.shape() ) ;
}

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value && boost::is_complex<typename A1::value_type>::value, unary_operation< ops_arg<typename A1::value_type>, A1 > >::type
arg( A1 const& a1 ) {
	return unary_operation< ops_arg<typename A1::value_type>, A1 >( ops_arg<typename A1::value_type>(), a1 )  ;
}

template <typename Y>
struct ops_cos {
	typedef decltype( std::cos( Y() ) ) result_type ;
	result_type operator() ( Y const& y ) const { return std::cos( y ) ; }
} ;

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value, unary_operation< ops_cos<typename A1::value_type>, A1 > >::type
cos( A1 const& a1 ) {
	return unary_operation< ops_cos<typename A1::value_type>, A1 >( ops_cos<typename A1::value_type>(), a1 )  ;
}

template <typename Y>
struct ops_acos {
	typedef decltype( std::acos( Y() ) ) result_type ;
	result_type operator() ( Y const& y ) const { return std::acos( y ) ; }
} ;

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value, unary_operation< ops_acos<typename A1::value_type>, A1 > >::type
acos( A1 const& a1 ) {
	return unary_operation< ops_acos<typename A1::value_type>, A1 >( ops_acos<typename A1::value_type>(), a1 )  ;
}

template <typename Y>
struct ops_cosh {
	typedef decltype( std::cosh( Y() ) ) result_type ;
	result_type operator() ( Y const& y ) const { return std::cosh( y ) ; }
} ;

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value, unary_operation< ops_cosh<typename A1::value_type>, A1 > >::type
cosh( A1 const& a1 ) {
	return unary_operation< ops_cosh<typename A1::value_type>, A1 >( ops_cosh<typename A1::value_type>(), a1 )  ;
}

template <typename Y>
struct ops_acosh {
	typedef decltype( std::acosh( Y() ) ) result_type ;
	result_type operator() ( Y const& y ) const { return std::acosh( y ) ; }
} ;

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value, unary_operation< ops_acosh<typename A1::value_type>, A1 > >::type
acosh( A1 const& a1 ) {
	return unary_operation< ops_acosh<typename A1::value_type>, A1 >( ops_acosh<typename A1::value_type>(), a1 )  ;
}

template <typename Y>
struct ops_sin {
	typedef decltype( std::sin( Y() ) ) result_type ;
	result_type operator() ( Y const& y ) const { return std::sin( y ) ; }
} ;

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value, unary_operation< ops_sin<typename A1::value_type>, A1 > >::type
sin( A1 const& a1 ) {
	return unary_operation< ops_sin<typename A1::value_type>, A1 >( ops_sin<typename A1::value_type>(), a1 )  ;
}

template <typename Y>
struct ops_asin {
	typedef decltype( std::asin( Y() ) ) result_type ;
	result_type operator() ( Y const& y ) const { return std::asin( y ) ; }
} ;

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value, unary_operation< ops_asin<typename A1::value_type>, A1 > >::type
asin( A1 const& a1 ) {
	return unary_operation< ops_asin<typename A1::value_type>, A1 >( ops_asin<typename A1::value_type>(), a1 )  ;
}

template <typename Y>
struct ops_sinh {
	typedef decltype( std::sinh( Y() ) ) result_type ;
	result_type operator() ( Y const& y ) const { return std::sinh( y ) ; }
} ;

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value, unary_operation< ops_sinh<typename A1::value_type>, A1 > >::type
sinh( A1 const& a1 ) {
	return unary_operation< ops_sinh<typename A1::value_type>, A1 >( ops_sinh<typename A1::value_type>(), a1 )  ;
}

template <typename Y>
struct ops_asinh {
	typedef decltype( std::asinh( Y() ) ) result_type ;
	result_type operator() ( Y const& y ) const { return std::asinh( y ) ; }
} ;

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value, unary_operation< ops_asinh<typename A1::value_type>, A1 > >::type
asinh( A1 const& a1 ) {
	return unary_operation< ops_asinh<typename A1::value_type>, A1 >( ops_asinh<typename A1::value_type>(), a1 )  ;
}

template <typename Y>
struct ops_tan {
	typedef decltype( std::tan( Y() ) ) result_type ;
	result_type operator() ( Y const& y ) const { return std::tan( y ) ; }
} ;

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value, unary_operation< ops_tan<typename A1::value_type>, A1 > >::type
tan( A1 const& a1 ) {
	return unary_operation< ops_tan<typename A1::value_type>, A1 >( ops_tan<typename A1::value_type>(), a1 )  ;
}

template <typename Y>
struct ops_atan {
	typedef decltype( std::atan( Y() ) ) result_type ;
	result_type operator() ( Y const& y ) const { return std::atan( y ) ; }
} ;

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value, unary_operation< ops_atan<typename A1::value_type>, A1 > >::type
atan( A1 const& a1 ) {
	return unary_operation< ops_atan<typename A1::value_type>, A1 >( ops_atan<typename A1::value_type>(), a1 )  ;
}

template <typename Y>
struct ops_tanh {
	typedef decltype( std::tanh( Y() ) ) result_type ;
	result_type operator() ( Y const& y ) const { return std::tanh( y ) ; }
} ;

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value, unary_operation< ops_tanh<typename A1::value_type>, A1 > >::type
tanh( A1 const& a1 ) {
	return unary_operation< ops_tanh<typename A1::value_type>, A1 >( ops_tanh<typename A1::value_type>(), a1 )  ;
}

template <typename Y>
struct ops_atanh {
	typedef decltype( std::atanh( Y() ) ) result_type ;
	result_type operator() ( Y const& y ) const { return std::atanh( y ) ; }
} ;

template <typename A1>
typename std::enable_if< is<DenseArray, A1>::value, unary_operation< ops_atanh<typename A1::value_type>, A1 > >::type
atanh( A1 const& a1 ) {
	return unary_operation< ops_atanh<typename A1::value_type>, A1 >( ops_atanh<typename A1::value_type>(), a1 )  ;
}

} // namespace glas3

#endif
