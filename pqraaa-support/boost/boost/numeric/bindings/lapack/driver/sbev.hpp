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

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_DRIVER_SBEV_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_DRIVER_SBEV_HPP

#include <cassert>
#include <boost/numeric/bindings/bandwidth.hpp>
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
// The LAPACK-backend for sbev is the netlib-compatible backend.
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
inline std::ptrdiff_t sbev( const char jobz, const UpLo uplo,
        const fortran_int_t n, const fortran_int_t kd, float* ab,
        const fortran_int_t ldab, float* w, float* z, const fortran_int_t ldz,
        float* work ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_SSBEV( &jobz, &lapack_option< UpLo >::value, &n, &kd, ab, &ldab, w,
            z, &ldz, work, &info );
#else
    LAPACK_SSBEV( &jobz, &lapack_option< UpLo >::value, &n, &kd, ab, &ldab, w,
            z, &ldz, work, &info ,1 ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * double value-type.
//
template< typename UpLo >
inline std::ptrdiff_t sbev( const char jobz, const UpLo uplo,
        const fortran_int_t n, const fortran_int_t kd, double* ab,
        const fortran_int_t ldab, double* w, double* z,
        const fortran_int_t ldz, double* work ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_DSBEV( &jobz, &lapack_option< UpLo >::value, &n, &kd, ab, &ldab, w,
            z, &ldz, work, &info );
#else
    LAPACK_DSBEV( &jobz, &lapack_option< UpLo >::value, &n, &kd, ab, &ldab, w,
            z, &ldz, work, &info ,1 ,1 );
#endif
    return info;
}

} // namespace detail

//
// Value-type based template class. Use this class if you need a type
// for dispatching to sbev.
//
template< typename Value >
struct sbev_impl {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function for user-defined workspaces, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename MatrixAB, typename VectorW, typename MatrixZ,
            typename WORK >
    static std::ptrdiff_t invoke( const char jobz, MatrixAB& ab, VectorW& w,
            MatrixZ& z, detail::workspace1< WORK > work ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::uplo_tag< MatrixAB >::type uplo;
        static_assert( (bindings::is_column_major< MatrixAB >::value) );
        static_assert( (bindings::is_column_major< MatrixZ >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixAB >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorW >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixAB >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixZ >::type >::type >::value) );
        static_assert( (bindings::is_mutable< MatrixAB >::value) );
        static_assert( (bindings::is_mutable< VectorW >::value) );
        static_assert( (bindings::is_mutable< MatrixZ >::value) );
        assert( bindings::bandwidth(ab, uplo()) >= 0 );
        assert( bindings::size(work.select(real_type())) >=
                min_size_work( bindings::size_column(ab) ));
        assert( bindings::size_column(ab) >= 0 );
        assert( bindings::size_minor(ab) == 1 ||
                bindings::stride_minor(ab) == 1 );
        assert( bindings::size_minor(z) == 1 ||
                bindings::stride_minor(z) == 1 );
        assert( bindings::stride_major(ab) >= bindings::bandwidth(ab,
                uplo())+1 );
        assert( jobz == 'N' || jobz == 'V' );
        return detail::sbev( jobz, uplo(), bindings::size_column(ab),
                bindings::bandwidth(ab, uplo()), bindings::begin_value(ab),
                bindings::stride_major(ab), bindings::begin_value(w),
                bindings::begin_value(z), bindings::stride_major(z),
                bindings::begin_value(work.select(real_type())) );
    }

    //
    // Static member function that
    // * Figures out the minimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member function
    // * Enables the unblocked algorithm (BLAS level 2)
    //
    template< typename MatrixAB, typename VectorW, typename MatrixZ >
    static std::ptrdiff_t invoke( const char jobz, MatrixAB& ab, VectorW& w,
            MatrixZ& z, minimal_workspace ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::uplo_tag< MatrixAB >::type uplo;
        bindings::detail::array< real_type > tmp_work( min_size_work(
                bindings::size_column(ab) ) );
        return invoke( jobz, ab, w, z, workspace( tmp_work ) );
    }

    //
    // Static member function that
    // * Figures out the optimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member
    // * Enables the blocked algorithm (BLAS level 3)
    //
    template< typename MatrixAB, typename VectorW, typename MatrixZ >
    static std::ptrdiff_t invoke( const char jobz, MatrixAB& ab, VectorW& w,
            MatrixZ& z, optimal_workspace ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::uplo_tag< MatrixAB >::type uplo;
        return invoke( jobz, ab, w, z, minimal_workspace() );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array work.
    //
    static std::ptrdiff_t min_size_work( const std::ptrdiff_t n ) {
        return std::max< std::ptrdiff_t >(1,3*n-2);
    }
};


//
// Functions for direct use. These functions are overloaded for temporaries,
// so that wrapped types can still be passed and used for write-access. In
// addition, if applicable, they are overloaded for user-defined workspaces.
// Calls to these functions are passed to the sbev_impl classes. In the 
// documentation, most overloads are collapsed to avoid a large number of
// prototypes which are very similar.
//

//
// Overloaded function for sbev. Its overload differs for
// * User-defined workspace
//
template< typename MatrixAB, typename VectorW, typename MatrixZ,
        typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
sbev( const char jobz, MatrixAB& ab, VectorW& w, MatrixZ& z,
        Workspace work ) {
    return sbev_impl< typename bindings::value_type<
            MatrixAB >::type >::invoke( jobz, ab, w, z, work );
}

//
// Overloaded function for sbev. Its overload differs for
// * Default workspace-type (optimal)
//
template< typename MatrixAB, typename VectorW, typename MatrixZ >
inline typename std::enable_if<! detail::is_workspace< MatrixZ >::value,
        std::ptrdiff_t >::type
sbev( const char jobz, MatrixAB& ab, VectorW& w, MatrixZ& z ) {
    return sbev_impl< typename bindings::value_type<
            MatrixAB >::type >::invoke( jobz, ab, w, z, optimal_workspace() );
}

} // namespace lapack
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
