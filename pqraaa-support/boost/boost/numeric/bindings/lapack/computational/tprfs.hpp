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

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_COMPUTATIONAL_TPRFS_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_COMPUTATIONAL_TPRFS_HPP

#include <cassert>
#include <boost/numeric/bindings/begin.hpp>
#include <boost/numeric/bindings/detail/array.hpp>
#include <boost/numeric/bindings/diag_tag.hpp>
#include <boost/numeric/bindings/is_column_major.hpp>
#include <boost/numeric/bindings/is_complex.hpp>
#include <boost/numeric/bindings/is_mutable.hpp>
#include <boost/numeric/bindings/is_real.hpp>
#include <boost/numeric/bindings/lapack/workspace.hpp>
#include <boost/numeric/bindings/remove_imaginary.hpp>
#include <boost/numeric/bindings/size.hpp>
#include <boost/numeric/bindings/stride.hpp>
#include <boost/numeric/bindings/uplo_tag.hpp>
#include <boost/numeric/bindings/value_type.hpp>
#include <type_traits>

//
// The LAPACK-backend for tprfs is the netlib-compatible backend.
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
inline std::ptrdiff_t tprfs( const UpLo uplo, const Trans trans,
        const Diag diag, const fortran_int_t n, const fortran_int_t nrhs,
        const float* ap, const float* b, const fortran_int_t ldb,
        const float* x, const fortran_int_t ldx, float* ferr, float* berr,
        float* work, fortran_int_t* iwork ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_STPRFS( &lapack_option< UpLo >::value, &lapack_option<
            Trans >::value, &lapack_option< Diag >::value, &n, &nrhs, ap, b,
            &ldb, x, &ldx, ferr, berr, work, iwork, &info );
#else
    LAPACK_STPRFS( &lapack_option< UpLo >::value, &lapack_option<
            Trans >::value, &lapack_option< Diag >::value, &n, &nrhs, ap, b,
            &ldb, x, &ldx, ferr, berr, work, iwork, &info ,1 ,1 ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * double value-type.
//
template< typename UpLo, typename Trans, typename Diag >
inline std::ptrdiff_t tprfs( const UpLo uplo, const Trans trans,
        const Diag diag, const fortran_int_t n, const fortran_int_t nrhs,
        const double* ap, const double* b, const fortran_int_t ldb,
        const double* x, const fortran_int_t ldx, double* ferr, double* berr,
        double* work, fortran_int_t* iwork ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_DTPRFS( &lapack_option< UpLo >::value, &lapack_option<
            Trans >::value, &lapack_option< Diag >::value, &n, &nrhs, ap, b,
            &ldb, x, &ldx, ferr, berr, work, iwork, &info );
#else
    LAPACK_DTPRFS( &lapack_option< UpLo >::value, &lapack_option<
            Trans >::value, &lapack_option< Diag >::value, &n, &nrhs, ap, b,
            &ldb, x, &ldx, ferr, berr, work, iwork, &info ,1 ,1 ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<float> value-type.
//
template< typename UpLo, typename Trans, typename Diag >
inline std::ptrdiff_t tprfs( const UpLo uplo, const Trans trans,
        const Diag diag, const fortran_int_t n, const fortran_int_t nrhs,
        const std::complex<float>* ap, const std::complex<float>* b,
        const fortran_int_t ldb, const std::complex<float>* x,
        const fortran_int_t ldx, float* ferr, float* berr,
        std::complex<float>* work, float* rwork ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_CTPRFS( &lapack_option< UpLo >::value, &lapack_option<
            Trans >::value, &lapack_option< Diag >::value, &n, &nrhs, ap, b,
            &ldb, x, &ldx, ferr, berr, work, rwork, &info );
#else
    LAPACK_CTPRFS( &lapack_option< UpLo >::value, &lapack_option<
            Trans >::value, &lapack_option< Diag >::value, &n, &nrhs, ap, b,
            &ldb, x, &ldx, ferr, berr, work, rwork, &info ,1 ,1 ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<double> value-type.
//
template< typename UpLo, typename Trans, typename Diag >
inline std::ptrdiff_t tprfs( const UpLo uplo, const Trans trans,
        const Diag diag, const fortran_int_t n, const fortran_int_t nrhs,
        const std::complex<double>* ap, const std::complex<double>* b,
        const fortran_int_t ldb, const std::complex<double>* x,
        const fortran_int_t ldx, double* ferr, double* berr,
        std::complex<double>* work, double* rwork ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_ZTPRFS( &lapack_option< UpLo >::value, &lapack_option<
            Trans >::value, &lapack_option< Diag >::value, &n, &nrhs, ap, b,
            &ldb, x, &ldx, ferr, berr, work, rwork, &info );
#else
    LAPACK_ZTPRFS( &lapack_option< UpLo >::value, &lapack_option<
            Trans >::value, &lapack_option< Diag >::value, &n, &nrhs, ap, b,
            &ldb, x, &ldx, ferr, berr, work, rwork, &info ,1 ,1 ,1 );
#endif
    return info;
}

} // namespace detail

//
// Value-type based template class. Use this class if you need a type
// for dispatching to tprfs.
//
template< typename Value, typename Enable = void >
struct tprfs_impl {};

//
// This implementation is enabled if Value is a real type.
//
template< typename Value >
struct tprfs_impl< Value, typename std::enable_if< is_real< Value >::value >::type > {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function for user-defined workspaces, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename MatrixAP, typename MatrixB, typename MatrixX,
            typename VectorFERR, typename VectorBERR, typename WORK,
            typename IWORK >
    static std::ptrdiff_t invoke( const MatrixAP& ap, const MatrixB& b,
            const MatrixX& x, VectorFERR& ferr, VectorBERR& berr,
            detail::workspace2< WORK, IWORK > work ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::uplo_tag< MatrixAP >::type uplo;
        typedef typename result_of::diag_tag< MatrixAP >::type diag;
        static_assert( (bindings::is_column_major< MatrixB >::value) );
        static_assert( (bindings::is_column_major< MatrixX >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixAP >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixB >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixAP >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixX >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixAP >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorFERR >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixAP >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorBERR >::type >::type >::value) );
        static_assert( (bindings::is_mutable< VectorFERR >::value) );
        static_assert( (bindings::is_mutable< VectorBERR >::value) );
        assert( bindings::size(berr) >= bindings::size_column(b) );
        assert( bindings::size(work.select(fortran_int_t())) >=
                min_size_iwork( bindings::size_column(ap) ));
        assert( bindings::size(work.select(real_type())) >=
                min_size_work( bindings::size_column(ap) ));
        assert( bindings::size_column(ap) >= 0 );
        assert( bindings::size_column(b) >= 0 );
        assert( bindings::size_minor(b) == 1 ||
                bindings::stride_minor(b) == 1 );
        assert( bindings::size_minor(x) == 1 ||
                bindings::stride_minor(x) == 1 );
        assert( bindings::stride_major(b) >= std::max< std::ptrdiff_t >(1,
                bindings::size_column(ap)) );
        assert( bindings::stride_major(x) >= std::max< std::ptrdiff_t >(1,
                bindings::size_column(ap)) );
        return detail::tprfs( uplo(), trans(), diag(),
                bindings::size_column(ap), bindings::size_column(b),
                bindings::begin_value(ap), bindings::begin_value(b),
                bindings::stride_major(b), bindings::begin_value(x),
                bindings::stride_major(x), bindings::begin_value(ferr),
                bindings::begin_value(berr),
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
    template< typename MatrixAP, typename MatrixB, typename MatrixX,
            typename VectorFERR, typename VectorBERR >
    static std::ptrdiff_t invoke( const MatrixAP& ap, const MatrixB& b,
            const MatrixX& x, VectorFERR& ferr, VectorBERR& berr,
            minimal_workspace work ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::uplo_tag< MatrixAP >::type uplo;
        typedef typename result_of::diag_tag< MatrixAP >::type diag;
        bindings::detail::array< real_type > tmp_work( min_size_work(
                bindings::size_column(ap) ) );
        bindings::detail::array< fortran_int_t > tmp_iwork(
                min_size_iwork( bindings::size_column(ap) ) );
        return invoke( ap, b, x, ferr, berr, workspace( tmp_work,
                tmp_iwork ) );
    }

    //
    // Static member function that
    // * Figures out the optimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member
    // * Enables the blocked algorithm (BLAS level 3)
    //
    template< typename MatrixAP, typename MatrixB, typename MatrixX,
            typename VectorFERR, typename VectorBERR >
    static std::ptrdiff_t invoke( const MatrixAP& ap, const MatrixB& b,
            const MatrixX& x, VectorFERR& ferr, VectorBERR& berr,
            optimal_workspace work ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::uplo_tag< MatrixAP >::type uplo;
        typedef typename result_of::diag_tag< MatrixAP >::type diag;
        return invoke( ap, b, x, ferr, berr, minimal_workspace() );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array work.
    //
    static std::ptrdiff_t min_size_work( const std::ptrdiff_t n ) {
        return 3*n;
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
struct tprfs_impl< Value, typename std::enable_if< is_complex< Value >::value >::type > {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function for user-defined workspaces, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename MatrixAP, typename MatrixB, typename MatrixX,
            typename VectorFERR, typename VectorBERR, typename WORK,
            typename RWORK >
    static std::ptrdiff_t invoke( const MatrixAP& ap, const MatrixB& b,
            const MatrixX& x, VectorFERR& ferr, VectorBERR& berr,
            detail::workspace2< WORK, RWORK > work ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::uplo_tag< MatrixAP >::type uplo;
        typedef typename result_of::diag_tag< MatrixAP >::type diag;
        static_assert( (bindings::is_column_major< MatrixB >::value) );
        static_assert( (bindings::is_column_major< MatrixX >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< VectorFERR >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorBERR >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixAP >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixB >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixAP >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixX >::type >::type >::value) );
        static_assert( (bindings::is_mutable< VectorFERR >::value) );
        static_assert( (bindings::is_mutable< VectorBERR >::value) );
        assert( bindings::size(berr) >= bindings::size_column(b) );
        assert( bindings::size(work.select(real_type())) >=
                min_size_rwork( bindings::size_column(ap) ));
        assert( bindings::size(work.select(value_type())) >=
                min_size_work( bindings::size_column(ap) ));
        assert( bindings::size_column(ap) >= 0 );
        assert( bindings::size_column(b) >= 0 );
        assert( bindings::size_minor(b) == 1 ||
                bindings::stride_minor(b) == 1 );
        assert( bindings::size_minor(x) == 1 ||
                bindings::stride_minor(x) == 1 );
        assert( bindings::stride_major(b) >= std::max< std::ptrdiff_t >(1,
                bindings::size_column(ap)) );
        assert( bindings::stride_major(x) >= std::max< std::ptrdiff_t >(1,
                bindings::size_column(ap)) );
        return detail::tprfs( uplo(), trans(), diag(),
                bindings::size_column(ap), bindings::size_column(b),
                bindings::begin_value(ap), bindings::begin_value(b),
                bindings::stride_major(b), bindings::begin_value(x),
                bindings::stride_major(x), bindings::begin_value(ferr),
                bindings::begin_value(berr),
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
    template< typename MatrixAP, typename MatrixB, typename MatrixX,
            typename VectorFERR, typename VectorBERR >
    static std::ptrdiff_t invoke( const MatrixAP& ap, const MatrixB& b,
            const MatrixX& x, VectorFERR& ferr, VectorBERR& berr,
            minimal_workspace work ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::uplo_tag< MatrixAP >::type uplo;
        typedef typename result_of::diag_tag< MatrixAP >::type diag;
        bindings::detail::array< value_type > tmp_work( min_size_work(
                bindings::size_column(ap) ) );
        bindings::detail::array< real_type > tmp_rwork( min_size_rwork(
                bindings::size_column(ap) ) );
        return invoke( ap, b, x, ferr, berr, workspace( tmp_work,
                tmp_rwork ) );
    }

    //
    // Static member function that
    // * Figures out the optimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member
    // * Enables the blocked algorithm (BLAS level 3)
    //
    template< typename MatrixAP, typename MatrixB, typename MatrixX,
            typename VectorFERR, typename VectorBERR >
    static std::ptrdiff_t invoke( const MatrixAP& ap, const MatrixB& b,
            const MatrixX& x, VectorFERR& ferr, VectorBERR& berr,
            optimal_workspace work ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::uplo_tag< MatrixAP >::type uplo;
        typedef typename result_of::diag_tag< MatrixAP >::type diag;
        return invoke( ap, b, x, ferr, berr, minimal_workspace() );
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
        return n;
    }
};


//
// Functions for direct use. These functions are overloaded for temporaries,
// so that wrapped types can still be passed and used for write-access. In
// addition, if applicable, they are overloaded for user-defined workspaces.
// Calls to these functions are passed to the tprfs_impl classes. In the 
// documentation, most overloads are collapsed to avoid a large number of
// prototypes which are very similar.
//

//
// Overloaded function for tprfs. Its overload differs for
// * VectorFERR&
// * VectorBERR&
// * User-defined workspace
//
template< typename MatrixAP, typename MatrixB, typename MatrixX,
        typename VectorFERR, typename VectorBERR, typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
tprfs( const MatrixAP& ap, const MatrixB& b, const MatrixX& x,
        VectorFERR& ferr, VectorBERR& berr, Workspace work ) {
    return tprfs_impl< typename bindings::value_type<
            MatrixAP >::type >::invoke( ap, b, x, ferr, berr, work );
}

//
// Overloaded function for tprfs. Its overload differs for
// * VectorFERR&
// * VectorBERR&
// * Default workspace-type (optimal)
//
template< typename MatrixAP, typename MatrixB, typename MatrixX,
        typename VectorFERR, typename VectorBERR >
inline typename std::enable_if<! detail::is_workspace< VectorBERR >::value,
        std::ptrdiff_t >::type
tprfs( const MatrixAP& ap, const MatrixB& b, const MatrixX& x,
        VectorFERR& ferr, VectorBERR& berr ) {
    return tprfs_impl< typename bindings::value_type<
            MatrixAP >::type >::invoke( ap, b, x, ferr, berr,
            optimal_workspace() );
}

//
// Overloaded function for tprfs. Its overload differs for
// * const VectorFERR&
// * VectorBERR&
// * User-defined workspace
//
template< typename MatrixAP, typename MatrixB, typename MatrixX,
        typename VectorFERR, typename VectorBERR, typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
tprfs( const MatrixAP& ap, const MatrixB& b, const MatrixX& x,
        const VectorFERR& ferr, VectorBERR& berr, Workspace work ) {
    return tprfs_impl< typename bindings::value_type<
            MatrixAP >::type >::invoke( ap, b, x, ferr, berr, work );
}

//
// Overloaded function for tprfs. Its overload differs for
// * const VectorFERR&
// * VectorBERR&
// * Default workspace-type (optimal)
//
template< typename MatrixAP, typename MatrixB, typename MatrixX,
        typename VectorFERR, typename VectorBERR >
inline typename std::enable_if<! detail::is_workspace< VectorBERR >::value,
        std::ptrdiff_t >::type
tprfs( const MatrixAP& ap, const MatrixB& b, const MatrixX& x,
        const VectorFERR& ferr, VectorBERR& berr ) {
    return tprfs_impl< typename bindings::value_type<
            MatrixAP >::type >::invoke( ap, b, x, ferr, berr,
            optimal_workspace() );
}

//
// Overloaded function for tprfs. Its overload differs for
// * VectorFERR&
// * const VectorBERR&
// * User-defined workspace
//
template< typename MatrixAP, typename MatrixB, typename MatrixX,
        typename VectorFERR, typename VectorBERR, typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
tprfs( const MatrixAP& ap, const MatrixB& b, const MatrixX& x,
        VectorFERR& ferr, const VectorBERR& berr, Workspace work ) {
    return tprfs_impl< typename bindings::value_type<
            MatrixAP >::type >::invoke( ap, b, x, ferr, berr, work );
}

//
// Overloaded function for tprfs. Its overload differs for
// * VectorFERR&
// * const VectorBERR&
// * Default workspace-type (optimal)
//
template< typename MatrixAP, typename MatrixB, typename MatrixX,
        typename VectorFERR, typename VectorBERR >
inline typename std::enable_if<! detail::is_workspace< VectorBERR >::value,
        std::ptrdiff_t >::type
tprfs( const MatrixAP& ap, const MatrixB& b, const MatrixX& x,
        VectorFERR& ferr, const VectorBERR& berr ) {
    return tprfs_impl< typename bindings::value_type<
            MatrixAP >::type >::invoke( ap, b, x, ferr, berr,
            optimal_workspace() );
}

//
// Overloaded function for tprfs. Its overload differs for
// * const VectorFERR&
// * const VectorBERR&
// * User-defined workspace
//
template< typename MatrixAP, typename MatrixB, typename MatrixX,
        typename VectorFERR, typename VectorBERR, typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
tprfs( const MatrixAP& ap, const MatrixB& b, const MatrixX& x,
        const VectorFERR& ferr, const VectorBERR& berr, Workspace work ) {
    return tprfs_impl< typename bindings::value_type<
            MatrixAP >::type >::invoke( ap, b, x, ferr, berr, work );
}

//
// Overloaded function for tprfs. Its overload differs for
// * const VectorFERR&
// * const VectorBERR&
// * Default workspace-type (optimal)
//
template< typename MatrixAP, typename MatrixB, typename MatrixX,
        typename VectorFERR, typename VectorBERR >
inline typename std::enable_if<! detail::is_workspace< VectorBERR >::value,
        std::ptrdiff_t >::type
tprfs( const MatrixAP& ap, const MatrixB& b, const MatrixX& x,
        const VectorFERR& ferr, const VectorBERR& berr ) {
    return tprfs_impl< typename bindings::value_type<
            MatrixAP >::type >::invoke( ap, b, x, ferr, berr,
            optimal_workspace() );
}

} // namespace lapack
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
