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

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_COMPUTATIONAL_TRTRI_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_COMPUTATIONAL_TRTRI_HPP

#include <cassert>
#include <boost/numeric/bindings/begin.hpp>
#include <boost/numeric/bindings/data_order.hpp>
#include <boost/numeric/bindings/diag_tag.hpp>
#include <boost/numeric/bindings/is_mutable.hpp>
#include <boost/numeric/bindings/remove_imaginary.hpp>
#include <boost/numeric/bindings/size.hpp>
#include <boost/numeric/bindings/stride.hpp>
#include <boost/numeric/bindings/uplo_tag.hpp>
#include <boost/numeric/bindings/value_type.hpp>
#include <type_traits>

//
// The LAPACK-backend for trtri is selected by defining a pre-processor
// variable, which can be one of
// * for ATLAS's CLAPACK, define BOOST_NUMERIC_BINDINGS_LAPACK_CLAPACK
// * netlib-compatible LAPACK is the default
//
#if defined BOOST_NUMERIC_BINDINGS_LAPACK_CLAPACK
#include <boost/numeric/bindings/lapack/detail/clapack.h>
#include <boost/numeric/bindings/lapack/detail/clapack_option.hpp>
#else
#include <boost/numeric/bindings/lapack/detail/lapack.h>
#include <boost/numeric/bindings/lapack/detail/lapack_option.hpp>
#endif

namespace boost {
namespace numeric {
namespace bindings {
namespace lapack {

//
// The detail namespace contains value-type-overloaded functions that
// dispatch to the appropriate back-end LAPACK-routine.
//
namespace detail {

#if defined BOOST_NUMERIC_BINDINGS_LAPACK_CLAPACK
//
// Overloaded function for dispatching to
// * ATLAS's CLAPACK backend, and
// * float value-type.
//
template< typename Order, typename UpLo, typename Diag >
inline std::ptrdiff_t trtri( Order, const UpLo uplo, const Diag diag,
        const int n, float* a, const int lda ) {
    return clapack_strtri( clapack_option< Order >::value, clapack_option<
            UpLo >::value, clapack_option< Diag >::value, n, a, lda );
}

//
// Overloaded function for dispatching to
// * ATLAS's CLAPACK backend, and
// * double value-type.
//
template< typename Order, typename UpLo, typename Diag >
inline std::ptrdiff_t trtri( Order, const UpLo uplo, const Diag diag,
        const int n, double* a, const int lda ) {
    return clapack_dtrtri( clapack_option< Order >::value, clapack_option<
            UpLo >::value, clapack_option< Diag >::value, n, a, lda );
}

//
// Overloaded function for dispatching to
// * ATLAS's CLAPACK backend, and
// * complex<float> value-type.
//
template< typename Order, typename UpLo, typename Diag >
inline std::ptrdiff_t trtri( Order, const UpLo uplo, const Diag diag,
        const int n, std::complex<float>* a, const int lda ) {
    return clapack_ctrtri( clapack_option< Order >::value, clapack_option<
            UpLo >::value, clapack_option< Diag >::value, n, a, lda );
}

//
// Overloaded function for dispatching to
// * ATLAS's CLAPACK backend, and
// * complex<double> value-type.
//
template< typename Order, typename UpLo, typename Diag >
inline std::ptrdiff_t trtri( Order, const UpLo uplo, const Diag diag,
        const int n, std::complex<double>* a, const int lda ) {
    return clapack_ztrtri( clapack_option< Order >::value, clapack_option<
            UpLo >::value, clapack_option< Diag >::value, n, a, lda );
}

#else
//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * float value-type.
//
template< typename Order, typename UpLo, typename Diag >
inline std::ptrdiff_t trtri( Order, const UpLo uplo, const Diag diag,
        const fortran_int_t n, float* a, const fortran_int_t lda ) {
    static_assert( (std::is_same<Order, tag::column_major>::value) );
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_STRTRI( &lapack_option< UpLo >::value, &lapack_option<
            Diag >::value, &n, a, &lda, &info );
#else
    LAPACK_STRTRI( &lapack_option< UpLo >::value, &lapack_option<
            Diag >::value, &n, a, &lda, &info ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * double value-type.
//
template< typename Order, typename UpLo, typename Diag >
inline std::ptrdiff_t trtri( Order, const UpLo uplo, const Diag diag,
        const fortran_int_t n, double* a, const fortran_int_t lda ) {
    static_assert( (std::is_same<Order, tag::column_major>::value) );
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_DTRTRI( &lapack_option< UpLo >::value, &lapack_option<
            Diag >::value, &n, a, &lda, &info );
#else
    LAPACK_DTRTRI( &lapack_option< UpLo >::value, &lapack_option<
            Diag >::value, &n, a, &lda, &info ,1 ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<float> value-type.
//
template< typename Order, typename UpLo, typename Diag >
inline std::ptrdiff_t trtri( Order, const UpLo uplo, const Diag diag,
        const fortran_int_t n, std::complex<float>* a,
        const fortran_int_t lda ) {
    static_assert( (std::is_same<Order, tag::column_major>::value) );
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_CTRTRI( &lapack_option< UpLo >::value, &lapack_option<
            Diag >::value, &n, a, &lda, &info );
#else
    LAPACK_CTRTRI( &lapack_option< UpLo >::value, &lapack_option<
            Diag >::value, &n, a, &lda, &info ,1 ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<double> value-type.
//
template< typename Order, typename UpLo, typename Diag >
inline std::ptrdiff_t trtri( Order, const UpLo uplo, const Diag diag,
        const fortran_int_t n, std::complex<double>* a,
        const fortran_int_t lda ) {
    static_assert( (std::is_same<Order, tag::column_major>::value) );
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_ZTRTRI( &lapack_option< UpLo >::value, &lapack_option<
            Diag >::value, &n, a, &lda, &info );
#else
    LAPACK_ZTRTRI( &lapack_option< UpLo >::value, &lapack_option<
            Diag >::value, &n, a, &lda, &info ,1 ,1 );
#endif
    return info;
}

#endif
} // namespace detail

//
// Value-type based template class. Use this class if you need a type
// for dispatching to trtri.
//
template< typename Value >
struct trtri_impl {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename MatrixA >
    static std::ptrdiff_t invoke( MatrixA& a ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::data_order< MatrixA >::type order;
        typedef typename result_of::uplo_tag< MatrixA >::type uplo;
        typedef typename result_of::diag_tag< MatrixA >::type diag;
        static_assert( (bindings::is_mutable< MatrixA >::value) );
        assert( bindings::size_column(a) >= 0 );
        assert( bindings::size_minor(a) == 1 ||
                bindings::stride_minor(a) == 1 );
        assert( bindings::stride_major(a) >= std::max< std::ptrdiff_t >(1,
                bindings::size_column(a)) );
        return detail::trtri( order(), uplo(), diag(),
                bindings::size_column(a), bindings::begin_value(a),
                bindings::stride_major(a) );
    }

};


//
// Functions for direct use. These functions are overloaded for temporaries,
// so that wrapped types can still be passed and used for write-access. In
// addition, if applicable, they are overloaded for user-defined workspaces.
// Calls to these functions are passed to the trtri_impl classes. In the 
// documentation, most overloads are collapsed to avoid a large number of
// prototypes which are very similar.
//

//
// Overloaded function for trtri. Its overload differs for
// * MatrixA&
//
template< typename MatrixA >
inline std::ptrdiff_t trtri( MatrixA& a ) {
    return trtri_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( a );
}

//
// Overloaded function for trtri. Its overload differs for
// * const MatrixA&
//
template< typename MatrixA >
inline std::ptrdiff_t trtri( const MatrixA& a ) {
    return trtri_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( a );
}

} // namespace lapack
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
