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

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_DRIVER_GGEVX_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_DRIVER_GGEVX_HPP

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
// The LAPACK-backend for ggevx is the netlib-compatible backend.
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
inline std::ptrdiff_t ggevx( const char balanc, const char jobvl,
        const char jobvr, const char sense, const fortran_int_t n, float* a,
        const fortran_int_t lda, float* b, const fortran_int_t ldb,
        float* alphar, float* alphai, float* beta, float* vl,
        const fortran_int_t ldvl, float* vr, const fortran_int_t ldvr,
        fortran_int_t& ilo, fortran_int_t& ihi, float* lscale, float* rscale,
        float& abnrm, float& bbnrm, float* rconde, float* rcondv, float* work,
        const fortran_int_t lwork, fortran_int_t* iwork,
        fortran_bool_t* bwork ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_SGGEVX( &balanc, &jobvl, &jobvr, &sense, &n, a, &lda, b, &ldb,
            alphar, alphai, beta, vl, &ldvl, vr, &ldvr, &ilo, &ihi, lscale,
            rscale, &abnrm, &bbnrm, rconde, rcondv, work, &lwork, iwork,
            bwork, &info );
#else
    LAPACK_SGGEVX( &balanc, &jobvl, &jobvr, &sense, &n, a, &lda, b, &ldb,
            alphar, alphai, beta, vl, &ldvl, vr, &ldvr, &ilo, &ihi, lscale,
            rscale, &abnrm, &bbnrm, rconde, rcondv, work, &lwork, iwork,
            bwork, &info ,1 ,1 ,1 ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * double value-type.
//
inline std::ptrdiff_t ggevx( const char balanc, const char jobvl,
        const char jobvr, const char sense, const fortran_int_t n, double* a,
        const fortran_int_t lda, double* b, const fortran_int_t ldb,
        double* alphar, double* alphai, double* beta, double* vl,
        const fortran_int_t ldvl, double* vr, const fortran_int_t ldvr,
        fortran_int_t& ilo, fortran_int_t& ihi, double* lscale,
        double* rscale, double& abnrm, double& bbnrm, double* rconde,
        double* rcondv, double* work, const fortran_int_t lwork,
        fortran_int_t* iwork, fortran_bool_t* bwork ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_DGGEVX( &balanc, &jobvl, &jobvr, &sense, &n, a, &lda, b, &ldb,
            alphar, alphai, beta, vl, &ldvl, vr, &ldvr, &ilo, &ihi, lscale,
            rscale, &abnrm, &bbnrm, rconde, rcondv, work, &lwork, iwork,
            bwork, &info );
#else
    LAPACK_DGGEVX( &balanc, &jobvl, &jobvr, &sense, &n, a, &lda, b, &ldb,
            alphar, alphai, beta, vl, &ldvl, vr, &ldvr, &ilo, &ihi, lscale,
            rscale, &abnrm, &bbnrm, rconde, rcondv, work, &lwork, iwork,
            bwork, &info ,1 ,1 ,1 ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<float> value-type.
//
inline std::ptrdiff_t ggevx( const char balanc, const char jobvl,
        const char jobvr, const char sense, const fortran_int_t n,
        std::complex<float>* a, const fortran_int_t lda,
        std::complex<float>* b, const fortran_int_t ldb,
        std::complex<float>* alpha, std::complex<float>* beta,
        std::complex<float>* vl, const fortran_int_t ldvl,
        std::complex<float>* vr, const fortran_int_t ldvr, fortran_int_t& ilo,
        fortran_int_t& ihi, float* lscale, float* rscale, float& abnrm,
        float& bbnrm, float* rconde, float* rcondv, std::complex<float>* work,
        const fortran_int_t lwork, float* rwork, fortran_int_t* iwork,
        fortran_bool_t* bwork ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_CGGEVX( &balanc, &jobvl, &jobvr, &sense, &n, a, &lda, b, &ldb,
            alpha, beta, vl, &ldvl, vr, &ldvr, &ilo, &ihi, lscale, rscale,
            &abnrm, &bbnrm, rconde, rcondv, work, &lwork, rwork, iwork, bwork,
            &info );
#else
    LAPACK_CGGEVX( &balanc, &jobvl, &jobvr, &sense, &n, a, &lda, b, &ldb,
            alpha, beta, vl, &ldvl, vr, &ldvr, &ilo, &ihi, lscale, rscale,
            &abnrm, &bbnrm, rconde, rcondv, work, &lwork, rwork, iwork, bwork,
            &info ,1 ,1 ,1 ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<double> value-type.
//
inline std::ptrdiff_t ggevx( const char balanc, const char jobvl,
        const char jobvr, const char sense, const fortran_int_t n,
        std::complex<double>* a, const fortran_int_t lda,
        std::complex<double>* b, const fortran_int_t ldb,
        std::complex<double>* alpha, std::complex<double>* beta,
        std::complex<double>* vl, const fortran_int_t ldvl,
        std::complex<double>* vr, const fortran_int_t ldvr,
        fortran_int_t& ilo, fortran_int_t& ihi, double* lscale,
        double* rscale, double& abnrm, double& bbnrm, double* rconde,
        double* rcondv, std::complex<double>* work, const fortran_int_t lwork,
        double* rwork, fortran_int_t* iwork, fortran_bool_t* bwork ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_ZGGEVX( &balanc, &jobvl, &jobvr, &sense, &n, a, &lda, b, &ldb,
            alpha, beta, vl, &ldvl, vr, &ldvr, &ilo, &ihi, lscale, rscale,
            &abnrm, &bbnrm, rconde, rcondv, work, &lwork, rwork, iwork, bwork,
            &info );
#else
    LAPACK_ZGGEVX( &balanc, &jobvl, &jobvr, &sense, &n, a, &lda, b, &ldb,
            alpha, beta, vl, &ldvl, vr, &ldvr, &ilo, &ihi, lscale, rscale,
            &abnrm, &bbnrm, rconde, rcondv, work, &lwork, rwork, iwork, bwork,
            &info ,1 ,1 ,1 ,1 );
#endif
    return info;
}

} // namespace detail

//
// Value-type based template class. Use this class if you need a type
// for dispatching to ggevx.
//
template< typename Value, typename Enable = void >
struct ggevx_impl {};

//
// This implementation is enabled if Value is a real type.
//
template< typename Value >
struct ggevx_impl< Value, typename std::enable_if< is_real< Value >::value >::type > {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function for user-defined workspaces, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename MatrixA, typename MatrixB, typename VectorALPHAR,
            typename VectorALPHAI, typename VectorBETA, typename MatrixVL,
            typename MatrixVR, typename VectorLSCALE, typename VectorRSCALE,
            typename VectorRCONDE, typename VectorRCONDV, typename WORK,
            typename IWORK, typename BWORK >
    static std::ptrdiff_t invoke( const char balanc, const char jobvl,
            const char jobvr, const char sense, MatrixA& a, MatrixB& b,
            VectorALPHAR& alphar, VectorALPHAI& alphai, VectorBETA& beta,
            MatrixVL& vl, MatrixVR& vr, fortran_int_t& ilo,
            fortran_int_t& ihi, VectorLSCALE& lscale,
            VectorRSCALE& rscale, real_type& abnrm, real_type& bbnrm,
            VectorRCONDE& rconde, VectorRCONDV& rcondv, detail::workspace3<
            WORK, IWORK, BWORK > work ) {
        namespace bindings = ::boost::numeric::bindings;
        static_assert( (bindings::is_column_major< MatrixA >::value) );
        static_assert( (bindings::is_column_major< MatrixB >::value) );
        static_assert( (bindings::is_column_major< MatrixVL >::value) );
        static_assert( (bindings::is_column_major< MatrixVR >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixB >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorALPHAR >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorALPHAI >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorBETA >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixVL >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixVR >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorLSCALE >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorRSCALE >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorRCONDE >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorRCONDV >::type >::type >::value) );
        static_assert( (bindings::is_mutable< MatrixA >::value) );
        static_assert( (bindings::is_mutable< MatrixB >::value) );
        static_assert( (bindings::is_mutable< VectorALPHAR >::value) );
        static_assert( (bindings::is_mutable< VectorALPHAI >::value) );
        static_assert( (bindings::is_mutable< VectorBETA >::value) );
        static_assert( (bindings::is_mutable< MatrixVL >::value) );
        static_assert( (bindings::is_mutable< MatrixVR >::value) );
        static_assert( (bindings::is_mutable< VectorLSCALE >::value) );
        static_assert( (bindings::is_mutable< VectorRSCALE >::value) );
        static_assert( (bindings::is_mutable< VectorRCONDE >::value) );
        static_assert( (bindings::is_mutable< VectorRCONDV >::value) );
        assert( bindings::size(alphai) >= bindings::size_column(a) );
        assert( bindings::size(alphar) >= bindings::size_column(a) );
        assert( bindings::size(work.select(fortran_int_t())) >=
                min_size_iwork( sense, bindings::size_column(a) ));
        assert( bindings::size(work.select(fortran_bool_t())) >=
                min_size_bwork( sense, bindings::size_column(a) ));
        assert( bindings::size(work.select(real_type())) >=
                min_size_work( balanc, jobvl, jobvr, sense,
                bindings::size_column(a) ));
        assert( bindings::size_column(a) >= 0 );
        assert( bindings::size_minor(a) == 1 ||
                bindings::stride_minor(a) == 1 );
        assert( bindings::size_minor(b) == 1 ||
                bindings::stride_minor(b) == 1 );
        assert( bindings::size_minor(vl) == 1 ||
                bindings::stride_minor(vl) == 1 );
        assert( bindings::size_minor(vr) == 1 ||
                bindings::stride_minor(vr) == 1 );
        assert( bindings::stride_major(a) >= std::max< std::ptrdiff_t >(1,
                bindings::size_column(a)) );
        assert( bindings::stride_major(b) >= std::max< std::ptrdiff_t >(1,
                bindings::size_column(a)) );
        assert( balanc == 'N' || balanc == 'P' || balanc == 'S' ||
                balanc == 'B' );
        assert( jobvl == 'N' || jobvl == 'V' );
        assert( jobvr == 'N' || jobvr == 'V' );
        assert( sense == 'N' || sense == 'E' || sense == 'V' ||
                sense == 'B' );
        return detail::ggevx( balanc, jobvl, jobvr, sense,
                bindings::size_column(a), bindings::begin_value(a),
                bindings::stride_major(a), bindings::begin_value(b),
                bindings::stride_major(b), bindings::begin_value(alphar),
                bindings::begin_value(alphai), bindings::begin_value(beta),
                bindings::begin_value(vl), bindings::stride_major(vl),
                bindings::begin_value(vr), bindings::stride_major(vr), ilo,
                ihi, bindings::begin_value(lscale),
                bindings::begin_value(rscale), abnrm, bbnrm,
                bindings::begin_value(rconde), bindings::begin_value(rcondv),
                bindings::begin_value(work.select(real_type())),
                bindings::size(work.select(real_type())),
                bindings::begin_value(work.select(fortran_int_t())),
                bindings::begin_value(work.select(fortran_bool_t())) );
    }

    //
    // Static member function that
    // * Figures out the minimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member function
    // * Enables the unblocked algorithm (BLAS level 2)
    //
    template< typename MatrixA, typename MatrixB, typename VectorALPHAR,
            typename VectorALPHAI, typename VectorBETA, typename MatrixVL,
            typename MatrixVR, typename VectorLSCALE, typename VectorRSCALE,
            typename VectorRCONDE, typename VectorRCONDV >
    static std::ptrdiff_t invoke( const char balanc, const char jobvl,
            const char jobvr, const char sense, MatrixA& a, MatrixB& b,
            VectorALPHAR& alphar, VectorALPHAI& alphai, VectorBETA& beta,
            MatrixVL& vl, MatrixVR& vr, fortran_int_t& ilo,
            fortran_int_t& ihi, VectorLSCALE& lscale,
            VectorRSCALE& rscale, real_type& abnrm, real_type& bbnrm,
            VectorRCONDE& rconde, VectorRCONDV& rcondv, minimal_workspace ) {
        namespace bindings = ::boost::numeric::bindings;
        bindings::detail::array< real_type > tmp_work( min_size_work( balanc,
                jobvl, jobvr, sense, bindings::size_column(a) ) );
        bindings::detail::array< fortran_int_t > tmp_iwork(
                min_size_iwork( sense, bindings::size_column(a) ) );
        bindings::detail::array< fortran_bool_t > tmp_bwork( min_size_bwork(
                sense, bindings::size_column(a) ) );
        return invoke( balanc, jobvl, jobvr, sense, a, b, alphar, alphai,
                beta, vl, vr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde,
                rcondv, workspace( tmp_work, tmp_iwork, tmp_bwork ) );
    }

    //
    // Static member function that
    // * Figures out the optimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member
    // * Enables the blocked algorithm (BLAS level 3)
    //
    template< typename MatrixA, typename MatrixB, typename VectorALPHAR,
            typename VectorALPHAI, typename VectorBETA, typename MatrixVL,
            typename MatrixVR, typename VectorLSCALE, typename VectorRSCALE,
            typename VectorRCONDE, typename VectorRCONDV >
    static std::ptrdiff_t invoke( const char balanc, const char jobvl,
            const char jobvr, const char sense, MatrixA& a, MatrixB& b,
            VectorALPHAR& alphar, VectorALPHAI& alphai, VectorBETA& beta,
            MatrixVL& vl, MatrixVR& vr, fortran_int_t& ilo,
            fortran_int_t& ihi, VectorLSCALE& lscale,
            VectorRSCALE& rscale, real_type& abnrm, real_type& bbnrm,
            VectorRCONDE& rconde, VectorRCONDV& rcondv, optimal_workspace ) {
        namespace bindings = ::boost::numeric::bindings;
        real_type opt_size_work;
        bindings::detail::array< fortran_int_t > tmp_iwork(
                min_size_iwork( sense, bindings::size_column(a) ) );
        bindings::detail::array< fortran_bool_t > tmp_bwork( min_size_bwork(
                sense, bindings::size_column(a) ) );
        detail::ggevx( balanc, jobvl, jobvr, sense,
                bindings::size_column(a), bindings::begin_value(a),
                bindings::stride_major(a), bindings::begin_value(b),
                bindings::stride_major(b), bindings::begin_value(alphar),
                bindings::begin_value(alphai), bindings::begin_value(beta),
                bindings::begin_value(vl), bindings::stride_major(vl),
                bindings::begin_value(vr), bindings::stride_major(vr), ilo,
                ihi, bindings::begin_value(lscale),
                bindings::begin_value(rscale), abnrm, bbnrm,
                bindings::begin_value(rconde), bindings::begin_value(rcondv),
                &opt_size_work, -1, bindings::begin_value(tmp_iwork),
                bindings::begin_value(tmp_bwork) );
        bindings::detail::array< real_type > tmp_work(
                traits::detail::to_int( opt_size_work ) );
        return invoke( balanc, jobvl, jobvr, sense, a, b, alphar, alphai,
                beta, vl, vr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde,
                rcondv, workspace( tmp_work, tmp_iwork, tmp_bwork ) );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array work.
    //
    static std::ptrdiff_t min_size_work( const char balanc, const char jobvl,
            const char jobvr, const char sense, const std::ptrdiff_t n ) {
        if ( balanc == 'S' || balanc == 'B' || jobvl == 'V' || jobvr == 'V' )
            return std::max< std::ptrdiff_t >( 1, 6*n );
        if ( sense == 'E' )
            return std::max< std::ptrdiff_t >( 1, 10*n );
        if ( sense == 'V' || sense == 'B' )
            return 2*n*n + 8*n + 16;
        return std::max< std::ptrdiff_t >( 1, 2*n );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array iwork.
    //
    static std::ptrdiff_t min_size_iwork( const char sense,
            const std::ptrdiff_t n ) {
        if ( sense == 'E' )
          return 0;
        else
          return n+6;
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array bwork.
    //
    static std::ptrdiff_t min_size_bwork( const char sense,
            const std::ptrdiff_t n ) {
        if ( sense == 'N' )
          return 0;
        else
          return n;
    }
};

//
// This implementation is enabled if Value is a complex type.
//
template< typename Value >
struct ggevx_impl< Value, typename std::enable_if< is_complex< Value >::value >::type > {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function for user-defined workspaces, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename MatrixA, typename MatrixB, typename VectorALPHA,
            typename VectorBETA, typename MatrixVL, typename MatrixVR,
            typename VectorLSCALE, typename VectorRSCALE,
            typename VectorRCONDE, typename VectorRCONDV, typename WORK,
            typename RWORK, typename IWORK, typename BWORK >
    static std::ptrdiff_t invoke( const char balanc, const char jobvl,
            const char jobvr, const char sense, MatrixA& a, MatrixB& b,
            VectorALPHA& alpha, VectorBETA& beta, MatrixVL& vl, MatrixVR& vr,
            fortran_int_t& ilo, fortran_int_t& ihi,
            VectorLSCALE& lscale, VectorRSCALE& rscale, real_type& abnrm,
            real_type& bbnrm, VectorRCONDE& rconde, VectorRCONDV& rcondv,
            detail::workspace4< WORK, RWORK, IWORK, BWORK > work ) {
        namespace bindings = ::boost::numeric::bindings;
        static_assert( (bindings::is_column_major< MatrixA >::value) );
        static_assert( (bindings::is_column_major< MatrixB >::value) );
        static_assert( (bindings::is_column_major< MatrixVL >::value) );
        static_assert( (bindings::is_column_major< MatrixVR >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< VectorLSCALE >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorRSCALE >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< VectorLSCALE >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorRCONDE >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< VectorLSCALE >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorRCONDV >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixB >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorALPHA >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorBETA >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixVL >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixVR >::type >::type >::value) );
        static_assert( (bindings::is_mutable< MatrixA >::value) );
        static_assert( (bindings::is_mutable< MatrixB >::value) );
        static_assert( (bindings::is_mutable< VectorALPHA >::value) );
        static_assert( (bindings::is_mutable< VectorBETA >::value) );
        static_assert( (bindings::is_mutable< MatrixVL >::value) );
        static_assert( (bindings::is_mutable< MatrixVR >::value) );
        static_assert( (bindings::is_mutable< VectorLSCALE >::value) );
        static_assert( (bindings::is_mutable< VectorRSCALE >::value) );
        static_assert( (bindings::is_mutable< VectorRCONDE >::value) );
        static_assert( (bindings::is_mutable< VectorRCONDV >::value) );
        assert( bindings::size(alpha) >= bindings::size_column(a) );
        assert( bindings::size(beta) >= bindings::size_column(a) );
        assert( bindings::size(work.select(fortran_int_t())) >=
                min_size_iwork( sense, bindings::size_column(a) ));
        assert( bindings::size(work.select(fortran_bool_t())) >=
                min_size_bwork( sense, bindings::size_column(a) ));
        assert( bindings::size(work.select(real_type())) >=
                min_size_rwork( balanc, bindings::size_column(a) ));
        assert( bindings::size(work.select(value_type())) >=
                min_size_work( sense, bindings::size_column(a) ));
        assert( bindings::size_column(a) >= 0 );
        assert( bindings::size_minor(a) == 1 ||
                bindings::stride_minor(a) == 1 );
        assert( bindings::size_minor(b) == 1 ||
                bindings::stride_minor(b) == 1 );
        assert( bindings::size_minor(vl) == 1 ||
                bindings::stride_minor(vl) == 1 );
        assert( bindings::size_minor(vr) == 1 ||
                bindings::stride_minor(vr) == 1 );
        assert( bindings::stride_major(a) >= std::max< std::ptrdiff_t >(1,
                bindings::size_column(a)) );
        assert( bindings::stride_major(b) >= std::max< std::ptrdiff_t >(1,
                bindings::size_column(a)) );
        assert( balanc == 'N' || balanc == 'P' || balanc == 'S' ||
                balanc == 'B' );
        assert( jobvl == 'N' || jobvl == 'V' );
        assert( jobvr == 'N' || jobvr == 'V' );
        assert( sense == 'N' || sense == 'E' || sense == 'V' ||
                sense == 'B' );
        return detail::ggevx( balanc, jobvl, jobvr, sense,
                bindings::size_column(a), bindings::begin_value(a),
                bindings::stride_major(a), bindings::begin_value(b),
                bindings::stride_major(b), bindings::begin_value(alpha),
                bindings::begin_value(beta), bindings::begin_value(vl),
                bindings::stride_major(vl), bindings::begin_value(vr),
                bindings::stride_major(vr), ilo, ihi,
                bindings::begin_value(lscale), bindings::begin_value(rscale),
                abnrm, bbnrm, bindings::begin_value(rconde),
                bindings::begin_value(rcondv),
                bindings::begin_value(work.select(value_type())),
                bindings::size(work.select(value_type())),
                bindings::begin_value(work.select(real_type())),
                bindings::begin_value(work.select(fortran_int_t())),
                bindings::begin_value(work.select(fortran_bool_t())) );
    }

    //
    // Static member function that
    // * Figures out the minimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member function
    // * Enables the unblocked algorithm (BLAS level 2)
    //
    template< typename MatrixA, typename MatrixB, typename VectorALPHA,
            typename VectorBETA, typename MatrixVL, typename MatrixVR,
            typename VectorLSCALE, typename VectorRSCALE,
            typename VectorRCONDE, typename VectorRCONDV >
    static std::ptrdiff_t invoke( const char balanc, const char jobvl,
            const char jobvr, const char sense, MatrixA& a, MatrixB& b,
            VectorALPHA& alpha, VectorBETA& beta, MatrixVL& vl, MatrixVR& vr,
            fortran_int_t& ilo, fortran_int_t& ihi,
            VectorLSCALE& lscale, VectorRSCALE& rscale, real_type& abnrm,
            real_type& bbnrm, VectorRCONDE& rconde, VectorRCONDV& rcondv,
            minimal_workspace ) {
        namespace bindings = ::boost::numeric::bindings;
        bindings::detail::array< value_type > tmp_work( min_size_work( sense,
                bindings::size_column(a) ) );
        bindings::detail::array< real_type > tmp_rwork( min_size_rwork(
                balanc, bindings::size_column(a) ) );
        bindings::detail::array< fortran_int_t > tmp_iwork(
                min_size_iwork( sense, bindings::size_column(a) ) );
        bindings::detail::array< fortran_bool_t > tmp_bwork( min_size_bwork(
                sense, bindings::size_column(a) ) );
        return invoke( balanc, jobvl, jobvr, sense, a, b, alpha, beta, vl, vr,
                ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv,
                workspace( tmp_work, tmp_rwork, tmp_iwork, tmp_bwork ) );
    }

    //
    // Static member function that
    // * Figures out the optimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member
    // * Enables the blocked algorithm (BLAS level 3)
    //
    template< typename MatrixA, typename MatrixB, typename VectorALPHA,
            typename VectorBETA, typename MatrixVL, typename MatrixVR,
            typename VectorLSCALE, typename VectorRSCALE,
            typename VectorRCONDE, typename VectorRCONDV >
    static std::ptrdiff_t invoke( const char balanc, const char jobvl,
            const char jobvr, const char sense, MatrixA& a, MatrixB& b,
            VectorALPHA& alpha, VectorBETA& beta, MatrixVL& vl, MatrixVR& vr,
            fortran_int_t& ilo, fortran_int_t& ihi,
            VectorLSCALE& lscale, VectorRSCALE& rscale, real_type& abnrm,
            real_type& bbnrm, VectorRCONDE& rconde, VectorRCONDV& rcondv,
            optimal_workspace ) {
        namespace bindings = ::boost::numeric::bindings;
        value_type opt_size_work;
        bindings::detail::array< real_type > tmp_rwork( min_size_rwork(
                balanc, bindings::size_column(a) ) );
        bindings::detail::array< fortran_int_t > tmp_iwork(
                min_size_iwork( sense, bindings::size_column(a) ) );
        bindings::detail::array< fortran_bool_t > tmp_bwork( min_size_bwork(
                sense, bindings::size_column(a) ) );
        detail::ggevx( balanc, jobvl, jobvr, sense,
                bindings::size_column(a), bindings::begin_value(a),
                bindings::stride_major(a), bindings::begin_value(b),
                bindings::stride_major(b), bindings::begin_value(alpha),
                bindings::begin_value(beta), bindings::begin_value(vl),
                bindings::stride_major(vl), bindings::begin_value(vr),
                bindings::stride_major(vr), ilo, ihi,
                bindings::begin_value(lscale), bindings::begin_value(rscale),
                abnrm, bbnrm, bindings::begin_value(rconde),
                bindings::begin_value(rcondv), &opt_size_work, -1,
                bindings::begin_value(tmp_rwork),
                bindings::begin_value(tmp_iwork),
                bindings::begin_value(tmp_bwork) );
        bindings::detail::array< value_type > tmp_work(
                traits::detail::to_int( opt_size_work ) );
        return invoke( balanc, jobvl, jobvr, sense, a, b, alpha, beta, vl, vr,
                ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv,
                workspace( tmp_work, tmp_rwork, tmp_iwork, tmp_bwork ) );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array work.
    //
    static std::ptrdiff_t min_size_work( const char sense,
            const std::ptrdiff_t n ) {
        if ( sense == 'N' )
            return std::max< std::ptrdiff_t >( 1, 2*n );
        else {
            if ( sense == 'E' )
                return std::max< std::ptrdiff_t >( 1, 4*n );
            else
                return std::max< std::ptrdiff_t >( 1, 2*n*n+2*n );
        }
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array rwork.
    //
    static std::ptrdiff_t min_size_rwork( const char balanc,
            const std::ptrdiff_t n ) {
        if ( balanc == 'S' || balanc == 'B' )
            return std::max< std::ptrdiff_t >( 1, 6*n );
        else
            return std::max< std::ptrdiff_t >( 1, 2*n );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array iwork.
    //
    static std::ptrdiff_t min_size_iwork( const char sense,
            const std::ptrdiff_t n ) {
        if ( sense == 'E' )
          return 0;
        else
          return n+2;
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array bwork.
    //
    static std::ptrdiff_t min_size_bwork( const char sense,
            const std::ptrdiff_t n ) {
        if ( sense == 'N' )
          return 0;
        else
          return n;
    }
};


//
// Functions for direct use. These functions are overloaded for temporaries,
// so that wrapped types can still be passed and used for write-access. In
// addition, if applicable, they are overloaded for user-defined workspaces.
// Calls to these functions are passed to the ggevx_impl classes. In the 
// documentation, most overloads are collapsed to avoid a large number of
// prototypes which are very similar.
//

//
// Overloaded function for ggevx. Its overload differs for
// * User-defined workspace
//
template< typename MatrixA, typename MatrixB, typename VectorALPHAR,
        typename VectorALPHAI, typename VectorBETA, typename MatrixVL,
        typename MatrixVR, typename VectorLSCALE, typename VectorRSCALE,
        typename VectorRCONDE, typename VectorRCONDV, typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
ggevx( const char balanc, const char jobvl, const char jobvr,
        const char sense, MatrixA& a, MatrixB& b, VectorALPHAR& alphar,
        VectorALPHAI& alphai, VectorBETA& beta, MatrixVL& vl, MatrixVR& vr,
        fortran_int_t& ilo, fortran_int_t& ihi, VectorLSCALE& lscale,
        VectorRSCALE& rscale, typename remove_imaginary<
        typename bindings::value_type< MatrixA >::type >::type& abnrm,
        typename remove_imaginary< typename bindings::value_type<
        MatrixA >::type >::type& bbnrm, VectorRCONDE& rconde,
        VectorRCONDV& rcondv, Workspace work ) {
    return ggevx_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( balanc, jobvl, jobvr, sense, a, b,
            alphar, alphai, beta, vl, vr, ilo, ihi, lscale, rscale, abnrm,
            bbnrm, rconde, rcondv, work );
}

//
// Overloaded function for ggevx. Its overload differs for
// * Default workspace-type (optimal)
//
template< typename MatrixA, typename MatrixB, typename VectorALPHAR,
        typename VectorALPHAI, typename VectorBETA, typename MatrixVL,
        typename MatrixVR, typename VectorLSCALE, typename VectorRSCALE,
        typename VectorRCONDE, typename VectorRCONDV >
inline typename std::enable_if<! detail::is_workspace< VectorRCONDV >::value,
        std::ptrdiff_t >::type
ggevx( const char balanc, const char jobvl, const char jobvr,
        const char sense, MatrixA& a, MatrixB& b, VectorALPHAR& alphar,
        VectorALPHAI& alphai, VectorBETA& beta, MatrixVL& vl, MatrixVR& vr,
        fortran_int_t& ilo, fortran_int_t& ihi, VectorLSCALE& lscale,
        VectorRSCALE& rscale, typename remove_imaginary<
        typename bindings::value_type< MatrixA >::type >::type& abnrm,
        typename remove_imaginary< typename bindings::value_type<
        MatrixA >::type >::type& bbnrm, VectorRCONDE& rconde,
        VectorRCONDV& rcondv ) {
    return ggevx_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( balanc, jobvl, jobvr, sense, a, b,
            alphar, alphai, beta, vl, vr, ilo, ihi, lscale, rscale, abnrm,
            bbnrm, rconde, rcondv, optimal_workspace() );
}
//
// Overloaded function for ggevx. Its overload differs for
// * User-defined workspace
//
template< typename MatrixA, typename MatrixB, typename VectorALPHA,
        typename VectorBETA, typename MatrixVL, typename MatrixVR,
        typename VectorLSCALE, typename VectorRSCALE, typename VectorRCONDE,
        typename VectorRCONDV, typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
ggevx( const char balanc, const char jobvl, const char jobvr,
        const char sense, MatrixA& a, MatrixB& b, VectorALPHA& alpha,
        VectorBETA& beta, MatrixVL& vl, MatrixVR& vr, fortran_int_t& ilo,
        fortran_int_t& ihi, VectorLSCALE& lscale, VectorRSCALE& rscale,
        typename remove_imaginary< typename bindings::value_type<
        MatrixA >::type >::type& abnrm, typename remove_imaginary<
        typename bindings::value_type< MatrixA >::type >::type& bbnrm,
        VectorRCONDE& rconde, VectorRCONDV& rcondv, Workspace work ) {
    return ggevx_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( balanc, jobvl, jobvr, sense, a, b,
            alpha, beta, vl, vr, ilo, ihi, lscale, rscale, abnrm, bbnrm,
            rconde, rcondv, work );
}

//
// Overloaded function for ggevx. Its overload differs for
// * Default workspace-type (optimal)
//
template< typename MatrixA, typename MatrixB, typename VectorALPHA,
        typename VectorBETA, typename MatrixVL, typename MatrixVR,
        typename VectorLSCALE, typename VectorRSCALE, typename VectorRCONDE,
        typename VectorRCONDV >
inline typename std::enable_if<! detail::is_workspace< VectorRCONDV >::value,
        std::ptrdiff_t >::type
ggevx( const char balanc, const char jobvl, const char jobvr,
        const char sense, MatrixA& a, MatrixB& b, VectorALPHA& alpha,
        VectorBETA& beta, MatrixVL& vl, MatrixVR& vr, fortran_int_t& ilo,
        fortran_int_t& ihi, VectorLSCALE& lscale, VectorRSCALE& rscale,
        typename remove_imaginary< typename bindings::value_type<
        MatrixA >::type >::type& abnrm, typename remove_imaginary<
        typename bindings::value_type< MatrixA >::type >::type& bbnrm,
        VectorRCONDE& rconde, VectorRCONDV& rcondv ) {
    return ggevx_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( balanc, jobvl, jobvr, sense, a, b,
            alpha, beta, vl, vr, ilo, ihi, lscale, rscale, abnrm, bbnrm,
            rconde, rcondv, optimal_workspace() );
}

} // namespace lapack
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
