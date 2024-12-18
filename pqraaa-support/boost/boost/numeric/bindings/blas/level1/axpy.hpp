//
// Copyright (c) 2002--2010
// Toon Knapen, Karl Meerbergen, Kresimir Fresl,
// Thomas Klimpel and Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
// THIS FILE IS AUTOMATICALLY GENERATED
// PLEASE DO NOT EDIT!
//

#ifndef BOOST_NUMERIC_BINDINGS_BLAS_LEVEL1_AXPY_HPP
#define BOOST_NUMERIC_BINDINGS_BLAS_LEVEL1_AXPY_HPP

#include <boost/numeric/bindings/begin.hpp>
#include <boost/numeric/bindings/is_mutable.hpp>
#include <boost/numeric/bindings/remove_imaginary.hpp>
#include <boost/numeric/bindings/size.hpp>
#include <boost/numeric/bindings/stride.hpp>
#include <boost/numeric/bindings/value_type.hpp>

//
// The BLAS-backend is selected by defining a pre-processor variable,
//  which can be one of
// * for CBLAS, define BOOST_NUMERIC_BINDINGS_BLAS_CBLAS
// * for CUBLAS, define BOOST_NUMERIC_BINDINGS_BLAS_CUBLAS
// * netlib-compatible BLAS is the default
//
#if defined BOOST_NUMERIC_BINDINGS_BLAS_CBLAS
#include <boost/numeric/bindings/blas/detail/cblas.h>
#include <boost/numeric/bindings/blas/detail/cblas_option.hpp>
#elif defined BOOST_NUMERIC_BINDINGS_BLAS_CUBLAS
#include <boost/numeric/bindings/blas/detail/cublas.h>
#include <boost/numeric/bindings/blas/detail/blas_option.hpp>
#else
#include <boost/numeric/bindings/blas/detail/blas.h>
#include <boost/numeric/bindings/blas/detail/blas_option.hpp>
#endif

namespace boost {
namespace numeric {
namespace bindings {
namespace blas {

//
// The detail namespace contains value-type-overloaded functions that
// dispatch to the appropriate back-end BLAS-routine.
//
namespace detail {

#if defined BOOST_NUMERIC_BINDINGS_BLAS_CBLAS
//
// Overloaded function for dispatching to
// * CBLAS backend, and
// * float value-type.
//
inline void axpy( const int n, const float a, const float* x, const int incx,
        float* y, const int incy ) {
    cblas_saxpy( n, a, x, incx, y, incy );
}

//
// Overloaded function for dispatching to
// * CBLAS backend, and
// * double value-type.
//
inline void axpy( const int n, const double a, const double* x,
        const int incx, double* y, const int incy ) {
    cblas_daxpy( n, a, x, incx, y, incy );
}

//
// Overloaded function for dispatching to
// * CBLAS backend, and
// * complex<float> value-type.
//
inline void axpy( const int n, const std::complex<float> a,
        const std::complex<float>* x, const int incx, std::complex<float>* y,
        const int incy ) {
    cblas_caxpy( n, &a, x, incx, y, incy );
}

//
// Overloaded function for dispatching to
// * CBLAS backend, and
// * complex<double> value-type.
//
inline void axpy( const int n, const std::complex<double> a,
        const std::complex<double>* x, const int incx,
        std::complex<double>* y, const int incy ) {
    cblas_zaxpy( n, &a, x, incx, y, incy );
}

#elif defined BOOST_NUMERIC_BINDINGS_BLAS_CUBLAS
//
// Overloaded function for dispatching to
// * CUBLAS backend, and
// * float value-type.
//
inline void axpy( const int n, const float a, const float* x, const int incx,
        float* y, const int incy ) {
    cublasSaxpy( n, a, x, incx, y, incy );
}

//
// Overloaded function for dispatching to
// * CUBLAS backend, and
// * double value-type.
//
inline void axpy( const int n, const double a, const double* x,
        const int incx, double* y, const int incy ) {
    cublasDaxpy( n, a, x, incx, y, incy );
}

//
// Overloaded function for dispatching to
// * CUBLAS backend, and
// * complex<float> value-type.
//
inline void axpy( const int n, const std::complex<float> a,
        const std::complex<float>* x, const int incx, std::complex<float>* y,
        const int incy ) {
    cublasCaxpy( n, a, x, incx, y, incy );
}

//
// Overloaded function for dispatching to
// * CUBLAS backend, and
// * complex<double> value-type.
//
inline void axpy( const int n, const std::complex<double> a,
        const std::complex<double>* x, const int incx,
        std::complex<double>* y, const int incy ) {
    // NOT FOUND();
}

#else
//
// Overloaded function for dispatching to
// * netlib-compatible BLAS backend (the default), and
// * float value-type.
//
inline void axpy( const fortran_int_t n, const float a, const float* x,
        const fortran_int_t incx, float* y, const fortran_int_t incy ) {
    BLAS_SAXPY( &n, &a, x, &incx, y, &incy );
}

//
// Overloaded function for dispatching to
// * netlib-compatible BLAS backend (the default), and
// * double value-type.
//
inline void axpy( const fortran_int_t n, const double a, const double* x,
        const fortran_int_t incx, double* y, const fortran_int_t incy ) {
    BLAS_DAXPY( &n, &a, x, &incx, y, &incy );
}

//
// Overloaded function for dispatching to
// * netlib-compatible BLAS backend (the default), and
// * complex<float> value-type.
//
inline void axpy( const fortran_int_t n, const std::complex<float> a,
        const std::complex<float>* x, const fortran_int_t incx,
        std::complex<float>* y, const fortran_int_t incy ) {
    BLAS_CAXPY( &n, &a, x, &incx, y, &incy );
}

//
// Overloaded function for dispatching to
// * netlib-compatible BLAS backend (the default), and
// * complex<double> value-type.
//
inline void axpy( const fortran_int_t n, const std::complex<double> a,
        const std::complex<double>* x, const fortran_int_t incx,
        std::complex<double>* y, const fortran_int_t incy ) {
    BLAS_ZAXPY( &n, &a, x, &incx, y, &incy );
}

#endif

} // namespace detail

//
// Value-type based template class. Use this class if you need a type
// for dispatching to axpy.
//
template< typename Value >
struct axpy_impl {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;
    typedef void return_type;

    //
    // Static member function that
    // * Deduces the required arguments for dispatching to BLAS, and
    // * Asserts that most arguments make sense.
    //
    template< typename VectorX, typename VectorY >
    static return_type invoke( const value_type a, const VectorX& x,
            VectorY& y ) {
        namespace bindings = ::boost::numeric::bindings;
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< VectorX >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorY >::type >::type >::value) );
        static_assert( (bindings::is_mutable< VectorY >::value) );
        detail::axpy( bindings::size(x), a, bindings::begin_value(x),
                bindings::stride(x), bindings::begin_value(y),
                bindings::stride(y) );
    }
};

//
// Functions for direct use. These functions are overloaded for temporaries,
// so that wrapped types can still be passed and used for write-access. Calls
// to these functions are passed to the axpy_impl classes. In the 
// documentation, the const-overloads are collapsed to avoid a large number of
// prototypes which are very similar.
//

//
// Overloaded function for axpy. Its overload differs for
// * VectorY&
//
template< typename VectorX, typename VectorY >
inline typename axpy_impl< typename bindings::value_type<
        VectorX >::type >::return_type
axpy( const typename bindings::value_type< VectorX >::type a,
        const VectorX& x, VectorY& y ) {
    axpy_impl< typename bindings::value_type<
            VectorX >::type >::invoke( a, x, y );
}

//
// Overloaded function for axpy. Its overload differs for
// * const VectorY&
//
template< typename VectorX, typename VectorY >
inline typename axpy_impl< typename bindings::value_type<
        VectorX >::type >::return_type
axpy( const typename bindings::value_type< VectorX >::type a,
        const VectorX& x, const VectorY& y ) {
    axpy_impl< typename bindings::value_type<
            VectorX >::type >::invoke( a, x, y );
}

} // namespace blas
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
