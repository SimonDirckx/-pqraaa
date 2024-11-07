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

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_DRIVER_SPEVD_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_DRIVER_SPEVD_HPP

#include <cassert>
#include <boost/numeric/bindings/begin.hpp>
#include <boost/numeric/bindings/detail/array.hpp>
#include <boost/numeric/bindings/is_column_major.hpp>
#include <boost/numeric/bindings/is_mutable.hpp>
#include <boost/numeric/bindings/lapack/workspace.hpp>
#include <boost/numeric/bindings/remove_imaginary.hpp>
#include <boost/numeric/bindings/size.hpp>
#include <boost/numeric/bindings/stride.hpp>
#include <boost/numeric/bindings/traits/detail/utils.hpp>
#include <boost/numeric/bindings/uplo_tag.hpp>
#include <boost/numeric/bindings/value_type.hpp>
#include <type_traits>

//
// The LAPACK-backend for spevd is the netlib-compatible backend.
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
inline std::ptrdiff_t spevd( const char jobz, const UpLo uplo,
        const fortran_int_t n, float* ap, float* w, float* z,
        const fortran_int_t ldz, float* work, const fortran_int_t lwork,
        fortran_int_t* iwork, const fortran_int_t liwork ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_SSPEVD( &jobz, &lapack_option< UpLo >::value, &n, ap, w, z, &ldz,
            work, &lwork, iwork, &liwork, &info );
#else
    LAPACK_SSPEVD( &jobz, &lapack_option< UpLo >::value, &n, ap, w, z, &ldz,
            work, &lwork, iwork, &liwork, &info ,1 ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * double value-type.
//
template< typename UpLo >
inline std::ptrdiff_t spevd( const char jobz, const UpLo uplo,
        const fortran_int_t n, double* ap, double* w, double* z,
        const fortran_int_t ldz, double* work, const fortran_int_t lwork,
        fortran_int_t* iwork, const fortran_int_t liwork ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_DSPEVD( &jobz, &lapack_option< UpLo >::value, &n, ap, w, z, &ldz,
            work, &lwork, iwork, &liwork, &info );
#else
    LAPACK_DSPEVD( &jobz, &lapack_option< UpLo >::value, &n, ap, w, z, &ldz,
            work, &lwork, iwork, &liwork, &info ,1 ,1 );
#endif
    return info;
}

} // namespace detail

//
// Value-type based template class. Use this class if you need a type
// for dispatching to spevd.
//
template< typename Value >
struct spevd_impl {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function for user-defined workspaces, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename MatrixAP, typename VectorW, typename MatrixZ,
            typename WORK, typename IWORK >
    static std::ptrdiff_t invoke( const char jobz, MatrixAP& ap, VectorW& w,
            MatrixZ& z, detail::workspace2< WORK, IWORK > work ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::uplo_tag< MatrixAP >::type uplo;
        static_assert( (bindings::is_column_major< MatrixZ >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixAP >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorW >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixAP >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixZ >::type >::type >::value) );
        static_assert( (bindings::is_mutable< MatrixAP >::value) );
        static_assert( (bindings::is_mutable< VectorW >::value) );
        static_assert( (bindings::is_mutable< MatrixZ >::value) );
        assert( bindings::size(work.select(fortran_int_t())) >=
                min_size_iwork( jobz, bindings::size_column(ap) ));
        assert( bindings::size(work.select(real_type())) >=
                min_size_work( jobz, bindings::size_column(ap) ));
        assert( bindings::size_column(ap) >= 0 );
        assert( bindings::size_minor(z) == 1 ||
                bindings::stride_minor(z) == 1 );
        assert( jobz == 'N' || jobz == 'V' );
        return detail::spevd( jobz, uplo(), bindings::size_column(ap),
                bindings::begin_value(ap), bindings::begin_value(w),
                bindings::begin_value(z), bindings::stride_major(z),
                bindings::begin_value(work.select(real_type())),
                bindings::size(work.select(real_type())),
                bindings::begin_value(work.select(fortran_int_t())),
                bindings::size(work.select(fortran_int_t())) );
    }

    //
    // Static member function that
    // * Figures out the minimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member function
    // * Enables the unblocked algorithm (BLAS level 2)
    //
    template< typename MatrixAP, typename VectorW, typename MatrixZ >
    static std::ptrdiff_t invoke( const char jobz, MatrixAP& ap, VectorW& w,
            MatrixZ& z, minimal_workspace ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::uplo_tag< MatrixAP >::type uplo;
        bindings::detail::array< real_type > tmp_work( min_size_work( jobz,
                bindings::size_column(ap) ) );
        bindings::detail::array< fortran_int_t > tmp_iwork(
                min_size_iwork( jobz, bindings::size_column(ap) ) );
        return invoke( jobz, ap, w, z, workspace( tmp_work, tmp_iwork ) );
    }

    //
    // Static member function that
    // * Figures out the optimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member
    // * Enables the blocked algorithm (BLAS level 3)
    //
    template< typename MatrixAP, typename VectorW, typename MatrixZ >
    static std::ptrdiff_t invoke( const char jobz, MatrixAP& ap, VectorW& w,
            MatrixZ& z, optimal_workspace ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::uplo_tag< MatrixAP >::type uplo;
        real_type opt_size_work;
        fortran_int_t opt_size_iwork;
        detail::spevd( jobz, uplo(), bindings::size_column(ap),
                bindings::begin_value(ap), bindings::begin_value(w),
                bindings::begin_value(z), bindings::stride_major(z),
                &opt_size_work, -1, &opt_size_iwork, -1 );
        bindings::detail::array< real_type > tmp_work(
                traits::detail::to_int( opt_size_work ) );
        bindings::detail::array< fortran_int_t > tmp_iwork(
                opt_size_iwork );
        return invoke( jobz, ap, w, z, workspace( tmp_work, tmp_iwork ) );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array work.
    //
    static std::ptrdiff_t min_size_work( const char jobz,
            const std::ptrdiff_t n ) {
        if ( n < 2 )
            return 1;
        else {
            if ( jobz == 'N' )
                return 2*n;
            else
                return 1 + 6*n + n*n;
        }
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array iwork.
    //
    static std::ptrdiff_t min_size_iwork( const char jobz,
            const std::ptrdiff_t n ) {
        if ( jobz == 'N' || n < 2 )
            return 1;
        else
            return 3 + 5*n;
    }
};


//
// Functions for direct use. These functions are overloaded for temporaries,
// so that wrapped types can still be passed and used for write-access. In
// addition, if applicable, they are overloaded for user-defined workspaces.
// Calls to these functions are passed to the spevd_impl classes. In the 
// documentation, most overloads are collapsed to avoid a large number of
// prototypes which are very similar.
//

//
// Overloaded function for spevd. Its overload differs for
// * User-defined workspace
//
template< typename MatrixAP, typename VectorW, typename MatrixZ,
        typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
spevd( const char jobz, MatrixAP& ap, VectorW& w, MatrixZ& z,
        Workspace work ) {
    return spevd_impl< typename bindings::value_type<
            MatrixAP >::type >::invoke( jobz, ap, w, z, work );
}

//
// Overloaded function for spevd. Its overload differs for
// * Default workspace-type (optimal)
//
template< typename MatrixAP, typename VectorW, typename MatrixZ >
inline typename std::enable_if<! detail::is_workspace< MatrixZ >::value,
        std::ptrdiff_t >::type
spevd( const char jobz, MatrixAP& ap, VectorW& w, MatrixZ& z ) {
    return spevd_impl< typename bindings::value_type<
            MatrixAP >::type >::invoke( jobz, ap, w, z, optimal_workspace() );
}

} // namespace lapack
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
