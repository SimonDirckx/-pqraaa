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

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_COMPUTATIONAL_TRTRS_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_COMPUTATIONAL_TRTRS_HPP

#include <cassert>
#include <boost/numeric/bindings/begin.hpp>
#include <boost/numeric/bindings/data_order.hpp>
#include <boost/numeric/bindings/diag_tag.hpp>
#include <boost/numeric/bindings/is_mutable.hpp>
#include <boost/numeric/bindings/remove_imaginary.hpp>
#include <boost/numeric/bindings/size.hpp>
#include <boost/numeric/bindings/stride.hpp>
#include <boost/numeric/bindings/trans_tag.hpp>
#include <boost/numeric/bindings/uplo_tag.hpp>
#include <boost/numeric/bindings/value_type.hpp>
#include <type_traits>

//
// The LAPACK-backend for trtrs is the netlib-compatible backend.
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
template< typename UpLo, typename Trans, typename Diag >
inline std::ptrdiff_t trtrs( const UpLo uplo, const Trans trans,
        const Diag diag, const fortran_int_t n, const fortran_int_t nrhs,
        const float* a, const fortran_int_t lda, float* b,
        const fortran_int_t ldb ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_STRTRS( &lapack_option< UpLo >::value, &lapack_option<
            Trans >::value, &lapack_option< Diag >::value, &n, &nrhs, a, &lda,
            b, &ldb, &info );
#else
    LAPACK_STRTRS( &lapack_option< UpLo >::value, &lapack_option<
            Trans >::value, &lapack_option< Diag >::value, &n, &nrhs, a, &lda,
            b, &ldb, &info ,1 ,1 ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * double value-type.
//
template< typename UpLo, typename Trans, typename Diag >
inline std::ptrdiff_t trtrs( const UpLo uplo, const Trans trans,
        const Diag diag, const fortran_int_t n, const fortran_int_t nrhs,
        const double* a, const fortran_int_t lda, double* b,
        const fortran_int_t ldb ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_DTRTRS( &lapack_option< UpLo >::value, &lapack_option<
            Trans >::value, &lapack_option< Diag >::value, &n, &nrhs, a, &lda,
            b, &ldb, &info );
#else
    LAPACK_DTRTRS( &lapack_option< UpLo >::value, &lapack_option<
            Trans >::value, &lapack_option< Diag >::value, &n, &nrhs, a, &lda,
            b, &ldb, &info ,1 ,1 ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<float> value-type.
//
template< typename UpLo, typename Trans, typename Diag >
inline std::ptrdiff_t trtrs( const UpLo uplo, const Trans trans,
        const Diag diag, const fortran_int_t n, const fortran_int_t nrhs,
        const std::complex<float>* a, const fortran_int_t lda,
        std::complex<float>* b, const fortran_int_t ldb ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_CTRTRS( &lapack_option< UpLo >::value, &lapack_option<
            Trans >::value, &lapack_option< Diag >::value, &n, &nrhs, a, &lda,
            b, &ldb, &info );
#else
    LAPACK_CTRTRS( &lapack_option< UpLo >::value, &lapack_option<
            Trans >::value, &lapack_option< Diag >::value, &n, &nrhs, a, &lda,
            b, &ldb, &info ,1 ,1 ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<double> value-type.
//
template< typename UpLo, typename Trans, typename Diag >
inline std::ptrdiff_t trtrs( const UpLo uplo, const Trans trans,
        const Diag diag, const fortran_int_t n, const fortran_int_t nrhs,
        const std::complex<double>* a, const fortran_int_t lda,
        std::complex<double>* b, const fortran_int_t ldb ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_ZTRTRS( &lapack_option< UpLo >::value, &lapack_option<
            Trans >::value, &lapack_option< Diag >::value, &n, &nrhs, a, &lda,
            b, &ldb, &info );
#else
    LAPACK_ZTRTRS( &lapack_option< UpLo >::value, &lapack_option<
            Trans >::value, &lapack_option< Diag >::value, &n, &nrhs, a, &lda,
            b, &ldb, &info ,1 ,1 ,1 );
#endif
    return info;
}

} // namespace detail

//
// Value-type based template class. Use this class if you need a type
// for dispatching to trtrs.
//
template< typename Value >
struct trtrs_impl {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename MatrixA, typename MatrixB >
    static std::ptrdiff_t invoke( const MatrixA& a, MatrixB& b ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::data_order< MatrixB >::type order;
        typedef typename result_of::trans_tag< MatrixA, order >::type trans;
        typedef typename result_of::uplo_tag< MatrixA, trans >::type uplo;
        typedef typename result_of::diag_tag< MatrixA >::type diag;
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixB >::type >::type >::value) );
        static_assert( (bindings::is_mutable< MatrixB >::value) );
        assert( bindings::size_column(b) >= 0 );
        assert( bindings::size_column_op(a, trans()) >= 0 );
        assert( bindings::size_minor(a) == 1 ||
                bindings::stride_minor(a) == 1 );
        assert( bindings::size_minor(b) == 1 ||
                bindings::stride_minor(b) == 1 );
        assert( bindings::stride_major(a) >= std::max< std::ptrdiff_t >(1,
                bindings::size_column_op(a, trans())) );
        assert( bindings::stride_major(b) >= std::max< std::ptrdiff_t >(1,
                bindings::size_column_op(a, trans())) );
        return detail::trtrs( uplo(), trans(), diag(),
                bindings::size_column_op(a, trans()),
                bindings::size_column(b), bindings::begin_value(a),
                bindings::stride_major(a), bindings::begin_value(b),
                bindings::stride_major(b) );
    }

};


//
// Functions for direct use. These functions are overloaded for temporaries,
// so that wrapped types can still be passed and used for write-access. In
// addition, if applicable, they are overloaded for user-defined workspaces.
// Calls to these functions are passed to the trtrs_impl classes. In the 
// documentation, most overloads are collapsed to avoid a large number of
// prototypes which are very similar.
//

//
// Overloaded function for trtrs. Its overload differs for
// * MatrixB&
//
template< typename MatrixA, typename MatrixB >
inline std::ptrdiff_t trtrs( const MatrixA& a, MatrixB& b ) {
    return trtrs_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( a, b );
}

//
// Overloaded function for trtrs. Its overload differs for
// * const MatrixB&
//
template< typename MatrixA, typename MatrixB >
inline std::ptrdiff_t trtrs( const MatrixA& a, const MatrixB& b ) {
    return trtrs_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( a, b );
}

} // namespace lapack
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
