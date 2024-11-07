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

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_DRIVER_GELSY_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_DRIVER_GELSY_HPP

#include <cassert>
#include <boost/numeric/bindings/begin.hpp>
#include <boost/numeric/bindings/detail/array.hpp>
#include <boost/numeric/bindings/is_column_major.hpp>
#include <boost/numeric/bindings/is_complex.hpp>
#include <boost/numeric/bindings/is_mutable.hpp>
#include <boost/numeric/bindings/is_real.hpp>
#include <boost/numeric/bindings/lapack/workspace.hpp>
#include <boost/numeric/bindings/remove_imaginary.hpp>
#include <boost/numeric/bindings/size.hpp>
#include <boost/numeric/bindings/stride.hpp>
#include <boost/numeric/bindings/traits/detail/utils.hpp>
#include <boost/numeric/bindings/value_type.hpp>
#include <type_traits>

//
// The LAPACK-backend for gelsy is the netlib-compatible backend.
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
inline std::ptrdiff_t gelsy( const fortran_int_t m, const fortran_int_t n,
        const fortran_int_t nrhs, float* a, const fortran_int_t lda, float* b,
        const fortran_int_t ldb, fortran_int_t* jpvt, const float rcond,
        fortran_int_t& rank, float* work, const fortran_int_t lwork ) {
    fortran_int_t info(0);
    LAPACK_SGELSY( &m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, &rank, work,
            &lwork, &info );
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * double value-type.
//
inline std::ptrdiff_t gelsy( const fortran_int_t m, const fortran_int_t n,
        const fortran_int_t nrhs, double* a, const fortran_int_t lda,
        double* b, const fortran_int_t ldb, fortran_int_t* jpvt,
        const double rcond, fortran_int_t& rank, double* work,
        const fortran_int_t lwork ) {
    fortran_int_t info(0);
    LAPACK_DGELSY( &m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, &rank, work,
            &lwork, &info );
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<float> value-type.
//
inline std::ptrdiff_t gelsy( const fortran_int_t m, const fortran_int_t n,
        const fortran_int_t nrhs, std::complex<float>* a,
        const fortran_int_t lda, std::complex<float>* b,
        const fortran_int_t ldb, fortran_int_t* jpvt, const float rcond,
        fortran_int_t& rank, std::complex<float>* work,
        const fortran_int_t lwork, float* rwork ) {
    fortran_int_t info(0);
    LAPACK_CGELSY( &m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, &rank, work,
            &lwork, rwork, &info );
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<double> value-type.
//
inline std::ptrdiff_t gelsy( const fortran_int_t m, const fortran_int_t n,
        const fortran_int_t nrhs, std::complex<double>* a,
        const fortran_int_t lda, std::complex<double>* b,
        const fortran_int_t ldb, fortran_int_t* jpvt, const double rcond,
        fortran_int_t& rank, std::complex<double>* work,
        const fortran_int_t lwork, double* rwork ) {
    fortran_int_t info(0);
    LAPACK_ZGELSY( &m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, &rank, work,
            &lwork, rwork, &info );
    return info;
}

} // namespace detail

//
// Value-type based template class. Use this class if you need a type
// for dispatching to gelsy.
//
template< typename Value, typename Enable = void >
struct gelsy_impl {};

//
// This implementation is enabled if Value is a real type.
//
template< typename Value >
struct gelsy_impl< Value, typename std::enable_if< is_real< Value >::value >::type > {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function for user-defined workspaces, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename MatrixA, typename MatrixB, typename VectorJPVT,
            typename WORK >
    static std::ptrdiff_t invoke( MatrixA& a, MatrixB& b, VectorJPVT& jpvt,
            const real_type rcond, fortran_int_t& rank,
            detail::workspace1< WORK > work ) {
        namespace bindings = ::boost::numeric::bindings;
        static_assert( (bindings::is_column_major< MatrixA >::value) );
        static_assert( (bindings::is_column_major< MatrixB >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixB >::type >::type >::value) );
        static_assert( (bindings::is_mutable< MatrixA >::value) );
        static_assert( (bindings::is_mutable< MatrixB >::value) );
        static_assert( (bindings::is_mutable< VectorJPVT >::value) );
        assert( bindings::size(work.select(real_type())) >=
                min_size_work( bindings::size_row(a),
                bindings::size_column(a), bindings::size_column(b) ));
        assert( bindings::size_column(a) >= 0 );
        assert( bindings::size_column(b) >= 0 );
        assert( bindings::size_minor(a) == 1 ||
                bindings::stride_minor(a) == 1 );
        assert( bindings::size_minor(b) == 1 ||
                bindings::stride_minor(b) == 1 );
        assert( bindings::size_row(a) >= 0 );
        assert( bindings::stride_major(a) >= std::max< std::ptrdiff_t >(1,
                bindings::size_row(a)) );
        assert( bindings::stride_major(b) >= std::max< std::ptrdiff_t >(1,
                std::max< std::ptrdiff_t >(bindings::size_row(a),
                bindings::size_column(a))) );
        return detail::gelsy( bindings::size_row(a), bindings::size_column(a),
                bindings::size_column(b), bindings::begin_value(a),
                bindings::stride_major(a), bindings::begin_value(b),
                bindings::stride_major(b), bindings::begin_value(jpvt), rcond,
                rank, bindings::begin_value(work.select(real_type())),
                bindings::size(work.select(real_type())) );
    }

    //
    // Static member function that
    // * Figures out the minimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member function
    // * Enables the unblocked algorithm (BLAS level 2)
    //
    template< typename MatrixA, typename MatrixB, typename VectorJPVT >
    static std::ptrdiff_t invoke( MatrixA& a, MatrixB& b, VectorJPVT& jpvt,
            const real_type rcond, fortran_int_t& rank,
            minimal_workspace ) {
        namespace bindings = ::boost::numeric::bindings;
        bindings::detail::array< real_type > tmp_work( min_size_work(
                bindings::size_row(a), bindings::size_column(a),
                bindings::size_column(b) ) );
        return invoke( a, b, jpvt, rcond, rank, workspace( tmp_work ) );
    }

    //
    // Static member function that
    // * Figures out the optimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member
    // * Enables the blocked algorithm (BLAS level 3)
    //
    template< typename MatrixA, typename MatrixB, typename VectorJPVT >
    static std::ptrdiff_t invoke( MatrixA& a, MatrixB& b, VectorJPVT& jpvt,
            const real_type rcond, fortran_int_t& rank,
            optimal_workspace ) {
        namespace bindings = ::boost::numeric::bindings;
        real_type opt_size_work;
        detail::gelsy( bindings::size_row(a), bindings::size_column(a),
                bindings::size_column(b), bindings::begin_value(a),
                bindings::stride_major(a), bindings::begin_value(b),
                bindings::stride_major(b), bindings::begin_value(jpvt), rcond,
                rank, &opt_size_work, -1 );
        bindings::detail::array< real_type > tmp_work(
                traits::detail::to_int( opt_size_work ) );
        return invoke( a, b, jpvt, rcond, rank, workspace( tmp_work ) );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array work.
    //
    static std::ptrdiff_t min_size_work( const std::ptrdiff_t m,
            const std::ptrdiff_t n, const std::ptrdiff_t nrhs ) {
        std::ptrdiff_t minmn = std::min< std::ptrdiff_t >( m, n );
        return std::max< std::ptrdiff_t >( 1, std::max< std::ptrdiff_t >( minmn+
                3*n+1, 2*minmn+nrhs ));
    }
};

//
// This implementation is enabled if Value is a complex type.
//
template< typename Value >
struct gelsy_impl< Value, typename std::enable_if< is_complex< Value >::value >::type > {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function for user-defined workspaces, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename MatrixA, typename MatrixB, typename VectorJPVT,
            typename WORK, typename RWORK >
    static std::ptrdiff_t invoke( MatrixA& a, MatrixB& b, VectorJPVT& jpvt,
            const real_type rcond, fortran_int_t& rank,
            detail::workspace2< WORK, RWORK > work ) {
        namespace bindings = ::boost::numeric::bindings;
        static_assert( (bindings::is_column_major< MatrixA >::value) );
        static_assert( (bindings::is_column_major< MatrixB >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixB >::type >::type >::value) );
        static_assert( (bindings::is_mutable< MatrixA >::value) );
        static_assert( (bindings::is_mutable< MatrixB >::value) );
        static_assert( (bindings::is_mutable< VectorJPVT >::value) );
        assert( bindings::size(work.select(real_type())) >=
                min_size_rwork( bindings::size_column(a) ));
        assert( bindings::size(work.select(value_type())) >=
                min_size_work( bindings::size_row(a),
                bindings::size_column(a), bindings::size_column(b) ));
        assert( bindings::size_column(a) >= 0 );
        assert( bindings::size_column(b) >= 0 );
        assert( bindings::size_minor(a) == 1 ||
                bindings::stride_minor(a) == 1 );
        assert( bindings::size_minor(b) == 1 ||
                bindings::stride_minor(b) == 1 );
        assert( bindings::size_row(a) >= 0 );
        assert( bindings::stride_major(a) >= std::max< std::ptrdiff_t >(1,
                bindings::size_row(a)) );
        assert( bindings::stride_major(b) >= std::max< std::ptrdiff_t >(1,
                std::max< std::ptrdiff_t >(bindings::size_row(a),
                bindings::size_column(a))) );
        return detail::gelsy( bindings::size_row(a), bindings::size_column(a),
                bindings::size_column(b), bindings::begin_value(a),
                bindings::stride_major(a), bindings::begin_value(b),
                bindings::stride_major(b), bindings::begin_value(jpvt), rcond,
                rank, bindings::begin_value(work.select(value_type())),
                bindings::size(work.select(value_type())),
                bindings::begin_value(work.select(real_type())) );
    }

    //
    // Static member function that
    // * Figures out the minimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member function
    // * Enables the unblocked algorithm (BLAS level 2)
    //
    template< typename MatrixA, typename MatrixB, typename VectorJPVT >
    static std::ptrdiff_t invoke( MatrixA& a, MatrixB& b, VectorJPVT& jpvt,
            const real_type rcond, fortran_int_t& rank,
            minimal_workspace ) {
        namespace bindings = ::boost::numeric::bindings;
        bindings::detail::array< value_type > tmp_work( min_size_work(
                bindings::size_row(a), bindings::size_column(a),
                bindings::size_column(b) ) );
        bindings::detail::array< real_type > tmp_rwork( min_size_rwork(
                bindings::size_column(a) ) );
        return invoke( a, b, jpvt, rcond, rank, workspace( tmp_work,
                tmp_rwork ) );
    }

    //
    // Static member function that
    // * Figures out the optimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member
    // * Enables the blocked algorithm (BLAS level 3)
    //
    template< typename MatrixA, typename MatrixB, typename VectorJPVT >
    static std::ptrdiff_t invoke( MatrixA& a, MatrixB& b, VectorJPVT& jpvt,
            const real_type rcond, fortran_int_t& rank,
            optimal_workspace ) {
        namespace bindings = ::boost::numeric::bindings;
        value_type opt_size_work;
        bindings::detail::array< real_type > tmp_rwork( min_size_rwork(
                bindings::size_column(a) ) );
        detail::gelsy( bindings::size_row(a), bindings::size_column(a),
                bindings::size_column(b), bindings::begin_value(a),
                bindings::stride_major(a), bindings::begin_value(b),
                bindings::stride_major(b), bindings::begin_value(jpvt), rcond,
                rank, &opt_size_work, -1, bindings::begin_value(tmp_rwork) );
        bindings::detail::array< value_type > tmp_work(
                traits::detail::to_int( opt_size_work ) );
        return invoke( a, b, jpvt, rcond, rank, workspace( tmp_work,
                tmp_rwork ) );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array work.
    //
    static std::ptrdiff_t min_size_work( const std::ptrdiff_t m,
            const std::ptrdiff_t n, const std::ptrdiff_t nrhs ) {
        std::ptrdiff_t minmn = std::min< std::ptrdiff_t >( m, n );
        return std::max< std::ptrdiff_t >( 1, std::max<
                std::ptrdiff_t >( std::max< std::ptrdiff_t >( 2*minmn, n+1 ),
                minmn+nrhs ) );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array rwork.
    //
    static std::ptrdiff_t min_size_rwork( const std::ptrdiff_t n ) {
        return 2*n;
    }
};


//
// Functions for direct use. These functions are overloaded for temporaries,
// so that wrapped types can still be passed and used for write-access. In
// addition, if applicable, they are overloaded for user-defined workspaces.
// Calls to these functions are passed to the gelsy_impl classes. In the 
// documentation, most overloads are collapsed to avoid a large number of
// prototypes which are very similar.
//

//
// Overloaded function for gelsy. Its overload differs for
// * User-defined workspace
//
template< typename MatrixA, typename MatrixB, typename VectorJPVT,
        typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
gelsy( MatrixA& a, MatrixB& b, VectorJPVT& jpvt,
        const typename remove_imaginary< typename bindings::value_type<
        MatrixA >::type >::type rcond, fortran_int_t& rank,
        Workspace work ) {
    return gelsy_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( a, b, jpvt, rcond, rank, work );
}

//
// Overloaded function for gelsy. Its overload differs for
// * Default workspace-type (optimal)
//
template< typename MatrixA, typename MatrixB, typename VectorJPVT >
inline typename std::enable_if<! detail::is_workspace< VectorJPVT >::value,
        std::ptrdiff_t >::type
gelsy( MatrixA& a, MatrixB& b, VectorJPVT& jpvt,
        const typename remove_imaginary< typename bindings::value_type<
        MatrixA >::type >::type rcond, fortran_int_t& rank ) {
    return gelsy_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( a, b, jpvt, rcond, rank,
            optimal_workspace() );
}

} // namespace lapack
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
