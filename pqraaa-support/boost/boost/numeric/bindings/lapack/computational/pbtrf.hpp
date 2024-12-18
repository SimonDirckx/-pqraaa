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

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_COMPUTATIONAL_PBTRF_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_COMPUTATIONAL_PBTRF_HPP

#include <cassert>
#include <boost/numeric/bindings/bandwidth.hpp>
#include <boost/numeric/bindings/begin.hpp>
#include <boost/numeric/bindings/is_column_major.hpp>
#include <boost/numeric/bindings/is_mutable.hpp>
#include <boost/numeric/bindings/remove_imaginary.hpp>
#include <boost/numeric/bindings/size.hpp>
#include <boost/numeric/bindings/stride.hpp>
#include <boost/numeric/bindings/uplo_tag.hpp>
#include <boost/numeric/bindings/value_type.hpp>
#include <type_traits>

//
// The LAPACK-backend for pbtrf is the netlib-compatible backend.
//
#include <boost/numeric/bindings/lapack/detail/lapack.h>
#include <boost/numeric/bindings/lapack/detail/lapack_option.hpp>

namespace boost {
namespace numeric {
namespace bindings {
namespace lapack {

//
// The detail namespace contains value-type-overloaded functions that
// dispatch to the appropriate back-end LAPACK-routine.
//
namespace detail {

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * float value-type.
//
template< typename UpLo >
inline std::ptrdiff_t pbtrf( const UpLo uplo, const fortran_int_t n,
        const fortran_int_t kd, float* ab, const fortran_int_t ldab ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_SPBTRF( &lapack_option< UpLo >::value, &n, &kd, ab, &ldab, &info );
#else
    LAPACK_SPBTRF( &lapack_option< UpLo >::value, &n, &kd, ab, &ldab, &info ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * double value-type.
//
template< typename UpLo >
inline std::ptrdiff_t pbtrf( const UpLo uplo, const fortran_int_t n,
        const fortran_int_t kd, double* ab, const fortran_int_t ldab ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_DPBTRF( &lapack_option< UpLo >::value, &n, &kd, ab, &ldab, &info );
#else
    LAPACK_DPBTRF( &lapack_option< UpLo >::value, &n, &kd, ab, &ldab, &info ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<float> value-type.
//
template< typename UpLo >
inline std::ptrdiff_t pbtrf( const UpLo uplo, const fortran_int_t n,
        const fortran_int_t kd, std::complex<float>* ab,
        const fortran_int_t ldab ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_CPBTRF( &lapack_option< UpLo >::value, &n, &kd, ab, &ldab, &info );
#else
    LAPACK_CPBTRF( &lapack_option< UpLo >::value, &n, &kd, ab, &ldab, &info ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<double> value-type.
//
template< typename UpLo >
inline std::ptrdiff_t pbtrf( const UpLo uplo, const fortran_int_t n,
        const fortran_int_t kd, std::complex<double>* ab,
        const fortran_int_t ldab ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_ZPBTRF( &lapack_option< UpLo >::value, &n, &kd, ab, &ldab, &info );
#else
    LAPACK_ZPBTRF( &lapack_option< UpLo >::value, &n, &kd, ab, &ldab, &info ,1 );
#endif
    return info;
}

} // namespace detail

//
// Value-type based template class. Use this class if you need a type
// for dispatching to pbtrf.
//
template< typename Value >
struct pbtrf_impl {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename MatrixAB >
    static std::ptrdiff_t invoke( MatrixAB& ab ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::uplo_tag< MatrixAB >::type uplo;
        static_assert( (bindings::is_column_major< MatrixAB >::value) );
        static_assert( (bindings::is_mutable< MatrixAB >::value) );
        assert( bindings::bandwidth(ab, uplo()) >= 0 );
        assert( bindings::size_column(ab) >= 0 );
        assert( bindings::size_minor(ab) == 1 ||
                bindings::stride_minor(ab) == 1 );
        assert( bindings::stride_major(ab) >= bindings::bandwidth(ab,
                uplo())+1 );
        return detail::pbtrf( uplo(), bindings::size_column(ab),
                bindings::bandwidth(ab, uplo()), bindings::begin_value(ab),
                bindings::stride_major(ab) );
    }

};


//
// Functions for direct use. These functions are overloaded for temporaries,
// so that wrapped types can still be passed and used for write-access. In
// addition, if applicable, they are overloaded for user-defined workspaces.
// Calls to these functions are passed to the pbtrf_impl classes. In the 
// documentation, most overloads are collapsed to avoid a large number of
// prototypes which are very similar.
//

//
// Overloaded function for pbtrf. Its overload differs for
// * MatrixAB&
//
template< typename MatrixAB >
inline std::ptrdiff_t pbtrf( MatrixAB& ab ) {
    return pbtrf_impl< typename bindings::value_type<
            MatrixAB >::type >::invoke( ab );
}

//
// Overloaded function for pbtrf. Its overload differs for
// * const MatrixAB&
//
template< typename MatrixAB >
inline std::ptrdiff_t pbtrf( const MatrixAB& ab ) {
    return pbtrf_impl< typename bindings::value_type<
            MatrixAB >::type >::invoke( ab );
}

} // namespace lapack
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
