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

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_DRIVER_SPGVX_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_DRIVER_SPGVX_HPP

#include <cassert>
#include <boost/numeric/bindings/begin.hpp>
#include <boost/numeric/bindings/detail/array.hpp>
#include <boost/numeric/bindings/is_column_major.hpp>
#include <boost/numeric/bindings/is_mutable.hpp>
#include <boost/numeric/bindings/lapack/workspace.hpp>
#include <boost/numeric/bindings/remove_imaginary.hpp>
#include <boost/numeric/bindings/size.hpp>
#include <boost/numeric/bindings/stride.hpp>
#include <boost/numeric/bindings/uplo_tag.hpp>
#include <boost/numeric/bindings/value_type.hpp>
#include <type_traits>

//
// The LAPACK-backend for spgvx is the netlib-compatible backend.
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
inline std::ptrdiff_t spgvx( const fortran_int_t itype, const char jobz,
        const char range, const UpLo uplo, const fortran_int_t n, float* ap,
        float* bp, const float vl, const float vu, const fortran_int_t il,
        const fortran_int_t iu, const float abstol, fortran_int_t& m,
        float* w, float* z, const fortran_int_t ldz, float* work,
        fortran_int_t* iwork, fortran_int_t* ifail ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_SSPGVX( &itype, &jobz, &range, &lapack_option< UpLo >::value, &n,
            ap, bp, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, work, iwork,
            ifail, &info );
#else
    LAPACK_SSPGVX( &itype, &jobz, &range, &lapack_option< UpLo >::value, &n,
            ap, bp, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, work, iwork,
            ifail, &info ,1 ,1 ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * double value-type.
//
template< typename UpLo >
inline std::ptrdiff_t spgvx( const fortran_int_t itype, const char jobz,
        const char range, const UpLo uplo, const fortran_int_t n, double* ap,
        double* bp, const double vl, const double vu, const fortran_int_t il,
        const fortran_int_t iu, const double abstol, fortran_int_t& m,
        double* w, double* z, const fortran_int_t ldz, double* work,
        fortran_int_t* iwork, fortran_int_t* ifail ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_DSPGVX( &itype, &jobz, &range, &lapack_option< UpLo >::value, &n,
            ap, bp, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, work, iwork,
            ifail, &info );
#else
    LAPACK_DSPGVX( &itype, &jobz, &range, &lapack_option< UpLo >::value, &n,
            ap, bp, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, work, iwork,
            ifail, &info ,1 ,1 ,1 );
#endif
    return info;
}

} // namespace detail

//
// Value-type based template class. Use this class if you need a type
// for dispatching to spgvx.
//
template< typename Value >
struct spgvx_impl {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function for user-defined workspaces, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename MatrixAP, typename MatrixBP, typename VectorW,
            typename MatrixZ, typename VectorIFAIL, typename WORK,
            typename IWORK >
    static std::ptrdiff_t invoke( const fortran_int_t itype,
            const char jobz, const char range, MatrixAP& ap, MatrixBP& bp,
            const real_type vl, const real_type vu,
            const fortran_int_t il, const fortran_int_t iu,
            const real_type abstol, fortran_int_t& m, VectorW& w,
            MatrixZ& z, VectorIFAIL& ifail, detail::workspace2< WORK,
            IWORK > work ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::uplo_tag< MatrixAP >::type uplo;
        static_assert( (bindings::is_column_major< MatrixZ >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixAP >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixBP >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixAP >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorW >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixAP >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixZ >::type >::type >::value) );
        static_assert( (bindings::is_mutable< MatrixAP >::value) );
        static_assert( (bindings::is_mutable< MatrixBP >::value) );
        static_assert( (bindings::is_mutable< VectorW >::value) );
        static_assert( (bindings::is_mutable< MatrixZ >::value) );
        static_assert( (bindings::is_mutable< VectorIFAIL >::value) );
        assert( bindings::size(w) >= bindings::size_column(ap) );
        assert( bindings::size(work.select(fortran_int_t())) >=
                min_size_iwork( bindings::size_column(ap) ));
        assert( bindings::size(work.select(real_type())) >=
                min_size_work( bindings::size_column(ap) ));
        assert( bindings::size_column(ap) >= 0 );
        assert( bindings::size_minor(z) == 1 ||
                bindings::stride_minor(z) == 1 );
        assert( jobz == 'N' || jobz == 'V' );
        assert( range == 'A' || range == 'V' || range == 'I' );
        return detail::spgvx( itype, jobz, range, uplo(),
                bindings::size_column(ap), bindings::begin_value(ap),
                bindings::begin_value(bp), vl, vu, il, iu, abstol, m,
                bindings::begin_value(w), bindings::begin_value(z),
                bindings::stride_major(z),
                bindings::begin_value(work.select(real_type())),
                bindings::begin_value(work.select(fortran_int_t())),
                bindings::begin_value(ifail) );
    }

    //
    // Static member function that
    // * Figures out the minimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member function
    // * Enables the unblocked algorithm (BLAS level 2)
    //
    template< typename MatrixAP, typename MatrixBP, typename VectorW,
            typename MatrixZ, typename VectorIFAIL >
    static std::ptrdiff_t invoke( const fortran_int_t itype,
            const char jobz, const char range, MatrixAP& ap, MatrixBP& bp,
            const real_type vl, const real_type vu,
            const fortran_int_t il, const fortran_int_t iu,
            const real_type abstol, fortran_int_t& m, VectorW& w,
            MatrixZ& z, VectorIFAIL& ifail, minimal_workspace ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::uplo_tag< MatrixAP >::type uplo;
        bindings::detail::array< real_type > tmp_work( min_size_work(
                bindings::size_column(ap) ) );
        bindings::detail::array< fortran_int_t > tmp_iwork(
                min_size_iwork( bindings::size_column(ap) ) );
        return invoke( itype, jobz, range, ap, bp, vl, vu, il, iu, abstol, m,
                w, z, ifail, workspace( tmp_work, tmp_iwork ) );
    }

    //
    // Static member function that
    // * Figures out the optimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member
    // * Enables the blocked algorithm (BLAS level 3)
    //
    template< typename MatrixAP, typename MatrixBP, typename VectorW,
            typename MatrixZ, typename VectorIFAIL >
    static std::ptrdiff_t invoke( const fortran_int_t itype,
            const char jobz, const char range, MatrixAP& ap, MatrixBP& bp,
            const real_type vl, const real_type vu,
            const fortran_int_t il, const fortran_int_t iu,
            const real_type abstol, fortran_int_t& m, VectorW& w,
            MatrixZ& z, VectorIFAIL& ifail, optimal_workspace ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::uplo_tag< MatrixAP >::type uplo;
        return invoke( itype, jobz, range, ap, bp, vl, vu, il, iu, abstol, m,
                w, z, ifail, minimal_workspace() );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array work.
    //
    static std::ptrdiff_t min_size_work( const std::ptrdiff_t n ) {
        return 8*n;
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array iwork.
    //
    static std::ptrdiff_t min_size_iwork( const std::ptrdiff_t n ) {
        return 5*n;
    }
};


//
// Functions for direct use. These functions are overloaded for temporaries,
// so that wrapped types can still be passed and used for write-access. In
// addition, if applicable, they are overloaded for user-defined workspaces.
// Calls to these functions are passed to the spgvx_impl classes. In the 
// documentation, most overloads are collapsed to avoid a large number of
// prototypes which are very similar.
//

//
// Overloaded function for spgvx. Its overload differs for
// * User-defined workspace
//
template< typename MatrixAP, typename MatrixBP, typename VectorW,
        typename MatrixZ, typename VectorIFAIL, typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
spgvx( const fortran_int_t itype, const char jobz, const char range,
        MatrixAP& ap, MatrixBP& bp, const typename remove_imaginary<
        typename bindings::value_type< MatrixAP >::type >::type vl,
        const typename remove_imaginary< typename bindings::value_type<
        MatrixAP >::type >::type vu, const fortran_int_t il,
        const fortran_int_t iu, const typename remove_imaginary<
        typename bindings::value_type< MatrixAP >::type >::type abstol,
        fortran_int_t& m, VectorW& w, MatrixZ& z, VectorIFAIL& ifail,
        Workspace work ) {
    return spgvx_impl< typename bindings::value_type<
            MatrixAP >::type >::invoke( itype, jobz, range, ap, bp, vl, vu,
            il, iu, abstol, m, w, z, ifail, work );
}

//
// Overloaded function for spgvx. Its overload differs for
// * Default workspace-type (optimal)
//
template< typename MatrixAP, typename MatrixBP, typename VectorW,
        typename MatrixZ, typename VectorIFAIL >
inline typename std::enable_if<! detail::is_workspace< VectorIFAIL >::value,
        std::ptrdiff_t >::type
spgvx( const fortran_int_t itype, const char jobz, const char range,
        MatrixAP& ap, MatrixBP& bp, const typename remove_imaginary<
        typename bindings::value_type< MatrixAP >::type >::type vl,
        const typename remove_imaginary< typename bindings::value_type<
        MatrixAP >::type >::type vu, const fortran_int_t il,
        const fortran_int_t iu, const typename remove_imaginary<
        typename bindings::value_type< MatrixAP >::type >::type abstol,
        fortran_int_t& m, VectorW& w, MatrixZ& z, VectorIFAIL& ifail ) {
    return spgvx_impl< typename bindings::value_type<
            MatrixAP >::type >::invoke( itype, jobz, range, ap, bp, vl, vu,
            il, iu, abstol, m, w, z, ifail, optimal_workspace() );
}

} // namespace lapack
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
