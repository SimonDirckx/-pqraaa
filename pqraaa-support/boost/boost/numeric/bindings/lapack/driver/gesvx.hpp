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

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_DRIVER_GESVX_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_DRIVER_GESVX_HPP

#include <cassert>
#include <boost/numeric/bindings/begin.hpp>
#include <boost/numeric/bindings/data_order.hpp>
#include <boost/numeric/bindings/detail/array.hpp>
#include <boost/numeric/bindings/is_column_major.hpp>
#include <boost/numeric/bindings/is_complex.hpp>
#include <boost/numeric/bindings/is_mutable.hpp>
#include <boost/numeric/bindings/is_real.hpp>
#include <boost/numeric/bindings/lapack/workspace.hpp>
#include <boost/numeric/bindings/remove_imaginary.hpp>
#include <boost/numeric/bindings/size.hpp>
#include <boost/numeric/bindings/stride.hpp>
#include <boost/numeric/bindings/trans_tag.hpp>
#include <boost/numeric/bindings/value_type.hpp>
#include <type_traits>

//
// The LAPACK-backend for gesvx is the netlib-compatible backend.
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
template< typename Trans >
inline std::ptrdiff_t gesvx( const char fact, const Trans trans,
        const fortran_int_t n, const fortran_int_t nrhs, float* a,
        const fortran_int_t lda, float* af, const fortran_int_t ldaf,
        fortran_int_t* ipiv, char& equed, float* r, float* c, float* b,
        const fortran_int_t ldb, float* x, const fortran_int_t ldx,
        float& rcond, float* ferr, float* berr, float* work,
        fortran_int_t* iwork ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_SGESVX( &fact, &lapack_option< Trans >::value, &n, &nrhs, a, &lda,
            af, &ldaf, ipiv, &equed, r, c, b, &ldb, x, &ldx, &rcond, ferr,
            berr, work, iwork, &info );
#else
    LAPACK_SGESVX( &fact, &lapack_option< Trans >::value, &n, &nrhs, a, &lda,
            af, &ldaf, ipiv, &equed, r, c, b, &ldb, x, &ldx, &rcond, ferr,
            berr, work, iwork, &info ,1 ,1 ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * double value-type.
//
template< typename Trans >
inline std::ptrdiff_t gesvx( const char fact, const Trans trans,
        const fortran_int_t n, const fortran_int_t nrhs, double* a,
        const fortran_int_t lda, double* af, const fortran_int_t ldaf,
        fortran_int_t* ipiv, char& equed, double* r, double* c, double* b,
        const fortran_int_t ldb, double* x, const fortran_int_t ldx,
        double& rcond, double* ferr, double* berr, double* work,
        fortran_int_t* iwork ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_DGESVX( &fact, &lapack_option< Trans >::value, &n, &nrhs, a, &lda,
            af, &ldaf, ipiv, &equed, r, c, b, &ldb, x, &ldx, &rcond, ferr,
            berr, work, iwork, &info );
#else
    LAPACK_DGESVX( &fact, &lapack_option< Trans >::value, &n, &nrhs, a, &lda,
            af, &ldaf, ipiv, &equed, r, c, b, &ldb, x, &ldx, &rcond, ferr,
            berr, work, iwork, &info ,1 ,1 ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<float> value-type.
//
template< typename Trans >
inline std::ptrdiff_t gesvx( const char fact, const Trans trans,
        const fortran_int_t n, const fortran_int_t nrhs,
        std::complex<float>* a, const fortran_int_t lda,
        std::complex<float>* af, const fortran_int_t ldaf,
        fortran_int_t* ipiv, char& equed, float* r, float* c,
        std::complex<float>* b, const fortran_int_t ldb,
        std::complex<float>* x, const fortran_int_t ldx, float& rcond,
        float* ferr, float* berr, std::complex<float>* work, float* rwork ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_CGESVX( &fact, &lapack_option< Trans >::value, &n, &nrhs, a, &lda,
            af, &ldaf, ipiv, &equed, r, c, b, &ldb, x, &ldx, &rcond, ferr,
            berr, work, rwork, &info );
#else
    LAPACK_CGESVX( &fact, &lapack_option< Trans >::value, &n, &nrhs, a, &lda,
            af, &ldaf, ipiv, &equed, r, c, b, &ldb, x, &ldx, &rcond, ferr,
            berr, work, rwork, &info ,1 ,1 ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<double> value-type.
//
template< typename Trans >
inline std::ptrdiff_t gesvx( const char fact, const Trans trans,
        const fortran_int_t n, const fortran_int_t nrhs,
        std::complex<double>* a, const fortran_int_t lda,
        std::complex<double>* af, const fortran_int_t ldaf,
        fortran_int_t* ipiv, char& equed, double* r, double* c,
        std::complex<double>* b, const fortran_int_t ldb,
        std::complex<double>* x, const fortran_int_t ldx, double& rcond,
        double* ferr, double* berr, std::complex<double>* work,
        double* rwork ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_ZGESVX( &fact, &lapack_option< Trans >::value, &n, &nrhs, a, &lda,
            af, &ldaf, ipiv, &equed, r, c, b, &ldb, x, &ldx, &rcond, ferr,
            berr, work, rwork, &info );
#else
    LAPACK_ZGESVX( &fact, &lapack_option< Trans >::value, &n, &nrhs, a, &lda,
            af, &ldaf, ipiv, &equed, r, c, b, &ldb, x, &ldx, &rcond, ferr,
            berr, work, rwork, &info ,1 ,1 ,1 );
#endif
    return info;
}

} // namespace detail

//
// Value-type based template class. Use this class if you need a type
// for dispatching to gesvx.
//
template< typename Value, typename Enable = void >
struct gesvx_impl {};

//
// This implementation is enabled if Value is a real type.
//
template< typename Value >
struct gesvx_impl< Value, typename std::enable_if< is_real< Value >::value >::type > {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function for user-defined workspaces, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename MatrixA, typename MatrixAF, typename VectorIPIV,
            typename VectorR, typename VectorC, typename MatrixB,
            typename MatrixX, typename VectorFERR, typename VectorBERR,
            typename WORK, typename IWORK >
    static std::ptrdiff_t invoke( const char fact, MatrixA& a, MatrixAF& af,
            VectorIPIV& ipiv, char& equed, VectorR& r, VectorC& c, MatrixB& b,
            MatrixX& x, real_type& rcond, VectorFERR& ferr, VectorBERR& berr,
            detail::workspace2< WORK, IWORK > work ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::data_order< MatrixAF >::type order;
        typedef typename result_of::trans_tag< MatrixA, order >::type trans;
        static_assert( (bindings::is_column_major< MatrixB >::value) );
        static_assert( (bindings::is_column_major< MatrixX >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixAF >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorR >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorC >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixB >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixX >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorFERR >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorBERR >::type >::type >::value) );
        static_assert( (bindings::is_mutable< MatrixA >::value) );
        static_assert( (bindings::is_mutable< MatrixAF >::value) );
        static_assert( (bindings::is_mutable< VectorIPIV >::value) );
        static_assert( (bindings::is_mutable< VectorR >::value) );
        static_assert( (bindings::is_mutable< VectorC >::value) );
        static_assert( (bindings::is_mutable< MatrixB >::value) );
        static_assert( (bindings::is_mutable< MatrixX >::value) );
        static_assert( (bindings::is_mutable< VectorFERR >::value) );
        static_assert( (bindings::is_mutable< VectorBERR >::value) );
        assert( bindings::size(berr) >= bindings::size_column(b) );
        assert( bindings::size(work.select(fortran_int_t())) >=
                min_size_iwork( bindings::size_column_op(a, trans()) ));
        assert( bindings::size(work.select(real_type())) >=
                min_size_work( bindings::size_column_op(a, trans()) ));
        assert( bindings::size_column(b) >= 0 );
        assert( bindings::size_column_op(a, trans()) >= 0 );
        assert( bindings::size_minor(a) == 1 ||
                bindings::stride_minor(a) == 1 );
        assert( bindings::size_minor(af) == 1 ||
                bindings::stride_minor(af) == 1 );
        assert( bindings::size_minor(b) == 1 ||
                bindings::stride_minor(b) == 1 );
        assert( bindings::size_minor(x) == 1 ||
                bindings::stride_minor(x) == 1 );
        assert( bindings::stride_major(a) >= std::max< std::ptrdiff_t >(1,
                bindings::size_column_op(a, trans())) );
        assert( bindings::stride_major(af) >= std::max<
                std::ptrdiff_t >(1,bindings::size_column_op(a, trans())) );
        assert( bindings::stride_major(b) >= std::max< std::ptrdiff_t >(1,
                bindings::size_column_op(a, trans())) );
        assert( bindings::stride_major(x) >= std::max< std::ptrdiff_t >(1,
                bindings::size_column_op(a, trans())) );
        assert( equed == 'N' || equed == 'R' || equed == 'C' ||
                equed == 'B' );
        assert( fact == 'F' || fact == 'N' || fact == 'E' );
        return detail::gesvx( fact, trans(), bindings::size_column_op(a,
                trans()), bindings::size_column(b), bindings::begin_value(a),
                bindings::stride_major(a), bindings::begin_value(af),
                bindings::stride_major(af), bindings::begin_value(ipiv),
                equed, bindings::begin_value(r), bindings::begin_value(c),
                bindings::begin_value(b), bindings::stride_major(b),
                bindings::begin_value(x), bindings::stride_major(x), rcond,
                bindings::begin_value(ferr), bindings::begin_value(berr),
                bindings::begin_value(work.select(real_type())),
                bindings::begin_value(work.select(fortran_int_t())) );
    }

    //
    // Static member function that
    // * Figures out the minimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member function
    // * Enables the unblocked algorithm (BLAS level 2)
    //
    template< typename MatrixA, typename MatrixAF, typename VectorIPIV,
            typename VectorR, typename VectorC, typename MatrixB,
            typename MatrixX, typename VectorFERR, typename VectorBERR >
    static std::ptrdiff_t invoke( const char fact, MatrixA& a, MatrixAF& af,
            VectorIPIV& ipiv, char& equed, VectorR& r, VectorC& c, MatrixB& b,
            MatrixX& x, real_type& rcond, VectorFERR& ferr, VectorBERR& berr,
            minimal_workspace ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::data_order< MatrixAF >::type order;
        typedef typename result_of::trans_tag< MatrixA, order >::type trans;
        bindings::detail::array< real_type > tmp_work( min_size_work(
                bindings::size_column_op(a, trans()) ) );
        bindings::detail::array< fortran_int_t > tmp_iwork(
                min_size_iwork( bindings::size_column_op(a, trans()) ) );
        return invoke( fact, a, af, ipiv, equed, r, c, b, x, rcond, ferr,
                berr, workspace( tmp_work, tmp_iwork ) );
    }

    //
    // Static member function that
    // * Figures out the optimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member
    // * Enables the blocked algorithm (BLAS level 3)
    //
    template< typename MatrixA, typename MatrixAF, typename VectorIPIV,
            typename VectorR, typename VectorC, typename MatrixB,
            typename MatrixX, typename VectorFERR, typename VectorBERR >
    static std::ptrdiff_t invoke( const char fact, MatrixA& a, MatrixAF& af,
            VectorIPIV& ipiv, char& equed, VectorR& r, VectorC& c, MatrixB& b,
            MatrixX& x, real_type& rcond, VectorFERR& ferr, VectorBERR& berr,
            optimal_workspace ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::data_order< MatrixAF >::type order;
        typedef typename result_of::trans_tag< MatrixA, order >::type trans;
        return invoke( fact, a, af, ipiv, equed, r, c, b, x, rcond, ferr,
                berr, minimal_workspace() );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array work.
    //
    static std::ptrdiff_t min_size_work( const std::ptrdiff_t n ) {
        return 4*n;
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array iwork.
    //
    static std::ptrdiff_t min_size_iwork( const std::ptrdiff_t n ) {
        return n;
    }
};

//
// This implementation is enabled if Value is a complex type.
//
template< typename Value >
struct gesvx_impl< Value, typename std::enable_if< is_complex< Value >::value >::type > {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function for user-defined workspaces, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename MatrixA, typename MatrixAF, typename VectorIPIV,
            typename VectorR, typename VectorC, typename MatrixB,
            typename MatrixX, typename VectorFERR, typename VectorBERR,
            typename WORK, typename RWORK >
    static std::ptrdiff_t invoke( const char fact, MatrixA& a, MatrixAF& af,
            VectorIPIV& ipiv, char& equed, VectorR& r, VectorC& c, MatrixB& b,
            MatrixX& x, real_type& rcond, VectorFERR& ferr, VectorBERR& berr,
            detail::workspace2< WORK, RWORK > work ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::data_order< MatrixAF >::type order;
        typedef typename result_of::trans_tag< MatrixA, order >::type trans;
        static_assert( (bindings::is_column_major< MatrixB >::value) );
        static_assert( (bindings::is_column_major< MatrixX >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< VectorR >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorC >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< VectorR >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorFERR >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< VectorR >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorBERR >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixAF >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixB >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixX >::type >::type >::value) );
        static_assert( (bindings::is_mutable< MatrixA >::value) );
        static_assert( (bindings::is_mutable< MatrixAF >::value) );
        static_assert( (bindings::is_mutable< VectorIPIV >::value) );
        static_assert( (bindings::is_mutable< VectorR >::value) );
        static_assert( (bindings::is_mutable< VectorC >::value) );
        static_assert( (bindings::is_mutable< MatrixB >::value) );
        static_assert( (bindings::is_mutable< MatrixX >::value) );
        static_assert( (bindings::is_mutable< VectorFERR >::value) );
        static_assert( (bindings::is_mutable< VectorBERR >::value) );
        assert( bindings::size(berr) >= bindings::size_column(b) );
        assert( bindings::size(work.select(real_type())) >=
                min_size_rwork( bindings::size_column_op(a, trans()) ));
        assert( bindings::size(work.select(value_type())) >=
                min_size_work( bindings::size_column_op(a, trans()) ));
        assert( bindings::size_column(b) >= 0 );
        assert( bindings::size_column_op(a, trans()) >= 0 );
        assert( bindings::size_minor(a) == 1 ||
                bindings::stride_minor(a) == 1 );
        assert( bindings::size_minor(af) == 1 ||
                bindings::stride_minor(af) == 1 );
        assert( bindings::size_minor(b) == 1 ||
                bindings::stride_minor(b) == 1 );
        assert( bindings::size_minor(x) == 1 ||
                bindings::stride_minor(x) == 1 );
        assert( bindings::stride_major(a) >= std::max< std::ptrdiff_t >(1,
                bindings::size_column_op(a, trans())) );
        assert( bindings::stride_major(af) >= std::max<
                std::ptrdiff_t >(1,bindings::size_column_op(a, trans())) );
        assert( bindings::stride_major(b) >= std::max< std::ptrdiff_t >(1,
                bindings::size_column_op(a, trans())) );
        assert( bindings::stride_major(x) >= std::max< std::ptrdiff_t >(1,
                bindings::size_column_op(a, trans())) );
        assert( equed == 'N' || equed == 'R' || equed == 'C' ||
                equed == 'B' );
        assert( fact == 'F' || fact == 'N' || fact == 'E' );
        return detail::gesvx( fact, trans(), bindings::size_column_op(a,
                trans()), bindings::size_column(b), bindings::begin_value(a),
                bindings::stride_major(a), bindings::begin_value(af),
                bindings::stride_major(af), bindings::begin_value(ipiv),
                equed, bindings::begin_value(r), bindings::begin_value(c),
                bindings::begin_value(b), bindings::stride_major(b),
                bindings::begin_value(x), bindings::stride_major(x), rcond,
                bindings::begin_value(ferr), bindings::begin_value(berr),
                bindings::begin_value(work.select(value_type())),
                bindings::begin_value(work.select(real_type())) );
    }

    //
    // Static member function that
    // * Figures out the minimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member function
    // * Enables the unblocked algorithm (BLAS level 2)
    //
    template< typename MatrixA, typename MatrixAF, typename VectorIPIV,
            typename VectorR, typename VectorC, typename MatrixB,
            typename MatrixX, typename VectorFERR, typename VectorBERR >
    static std::ptrdiff_t invoke( const char fact, MatrixA& a, MatrixAF& af,
            VectorIPIV& ipiv, char& equed, VectorR& r, VectorC& c, MatrixB& b,
            MatrixX& x, real_type& rcond, VectorFERR& ferr, VectorBERR& berr,
            minimal_workspace ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::data_order< MatrixAF >::type order;
        typedef typename result_of::trans_tag< MatrixA, order >::type trans;
        bindings::detail::array< value_type > tmp_work( min_size_work(
                bindings::size_column_op(a, trans()) ) );
        bindings::detail::array< real_type > tmp_rwork( min_size_rwork(
                bindings::size_column_op(a, trans()) ) );
        return invoke( fact, a, af, ipiv, equed, r, c, b, x, rcond, ferr,
                berr, workspace( tmp_work, tmp_rwork ) );
    }

    //
    // Static member function that
    // * Figures out the optimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member
    // * Enables the blocked algorithm (BLAS level 3)
    //
    template< typename MatrixA, typename MatrixAF, typename VectorIPIV,
            typename VectorR, typename VectorC, typename MatrixB,
            typename MatrixX, typename VectorFERR, typename VectorBERR >
    static std::ptrdiff_t invoke( const char fact, MatrixA& a, MatrixAF& af,
            VectorIPIV& ipiv, char& equed, VectorR& r, VectorC& c, MatrixB& b,
            MatrixX& x, real_type& rcond, VectorFERR& ferr, VectorBERR& berr,
            optimal_workspace ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::data_order< MatrixAF >::type order;
        typedef typename result_of::trans_tag< MatrixA, order >::type trans;
        return invoke( fact, a, af, ipiv, equed, r, c, b, x, rcond, ferr,
                berr, minimal_workspace() );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array work.
    //
    static std::ptrdiff_t min_size_work( const std::ptrdiff_t n ) {
        return 2*n;
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
// Calls to these functions are passed to the gesvx_impl classes. In the 
// documentation, most overloads are collapsed to avoid a large number of
// prototypes which are very similar.
//

//
// Overloaded function for gesvx. Its overload differs for
// * User-defined workspace
//
template< typename MatrixA, typename MatrixAF, typename VectorIPIV,
        typename VectorR, typename VectorC, typename MatrixB,
        typename MatrixX, typename VectorFERR, typename VectorBERR,
        typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
gesvx( const char fact, MatrixA& a, MatrixAF& af, VectorIPIV& ipiv,
        char& equed, VectorR& r, VectorC& c, MatrixB& b, MatrixX& x,
        typename remove_imaginary< typename bindings::value_type<
        MatrixA >::type >::type& rcond, VectorFERR& ferr, VectorBERR& berr,
        Workspace work ) {
    return gesvx_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( fact, a, af, ipiv, equed, r, c, b, x,
            rcond, ferr, berr, work );
}

//
// Overloaded function for gesvx. Its overload differs for
// * Default workspace-type (optimal)
//
template< typename MatrixA, typename MatrixAF, typename VectorIPIV,
        typename VectorR, typename VectorC, typename MatrixB,
        typename MatrixX, typename VectorFERR, typename VectorBERR >
inline typename std::enable_if<! detail::is_workspace< VectorBERR >::value,
        std::ptrdiff_t >::type
gesvx( const char fact, MatrixA& a, MatrixAF& af, VectorIPIV& ipiv,
        char& equed, VectorR& r, VectorC& c, MatrixB& b, MatrixX& x,
        typename remove_imaginary< typename bindings::value_type<
        MatrixA >::type >::type& rcond, VectorFERR& ferr, VectorBERR& berr ) {
    return gesvx_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( fact, a, af, ipiv, equed, r, c, b, x,
            rcond, ferr, berr, optimal_workspace() );
}

} // namespace lapack
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
