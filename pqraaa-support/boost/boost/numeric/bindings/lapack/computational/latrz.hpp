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

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_COMPUTATIONAL_LATRZ_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_COMPUTATIONAL_LATRZ_HPP

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
#include <boost/numeric/bindings/value_type.hpp>
#include <type_traits>

//
// The LAPACK-backend for latrz is the netlib-compatible backend.
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
inline std::ptrdiff_t latrz( const fortran_int_t m, const fortran_int_t n,
        const fortran_int_t l, float* a, const fortran_int_t lda, float* tau,
        float* work ) {
    fortran_int_t info(0);
    LAPACK_SLATRZ( &m, &n, &l, a, &lda, tau, work );
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * double value-type.
//
inline std::ptrdiff_t latrz( const fortran_int_t m, const fortran_int_t n,
        const fortran_int_t l, double* a, const fortran_int_t lda,
        double* tau, double* work ) {
    fortran_int_t info(0);
    LAPACK_DLATRZ( &m, &n, &l, a, &lda, tau, work );
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<float> value-type.
//
inline std::ptrdiff_t latrz( const fortran_int_t m, const fortran_int_t n,
        const fortran_int_t l, std::complex<float>* a,
        const fortran_int_t lda, std::complex<float>* tau,
        std::complex<float>* work ) {
    fortran_int_t info(0);
    LAPACK_CLATRZ( &m, &n, &l, a, &lda, tau, work );
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<double> value-type.
//
inline std::ptrdiff_t latrz( const fortran_int_t m, const fortran_int_t n,
        const fortran_int_t l, std::complex<double>* a,
        const fortran_int_t lda, std::complex<double>* tau,
        std::complex<double>* work ) {
    fortran_int_t info(0);
    LAPACK_ZLATRZ( &m, &n, &l, a, &lda, tau, work );
    return info;
}

} // namespace detail

//
// Value-type based template class. Use this class if you need a type
// for dispatching to latrz.
//
template< typename Value, typename Enable = void >
struct latrz_impl {};

//
// This implementation is enabled if Value is a real type.
//
template< typename Value >
struct latrz_impl< Value, typename std::enable_if< is_real< Value >::value >::type > {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function for user-defined workspaces, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename MatrixA, typename VectorTAU, typename WORK >
    static std::ptrdiff_t invoke( MatrixA& a, VectorTAU& tau,
            detail::workspace1< WORK > work ) {
        namespace bindings = ::boost::numeric::bindings;
        static_assert( (bindings::is_column_major< MatrixA >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorTAU >::type >::type >::value) );
        static_assert( (bindings::is_mutable< MatrixA >::value) );
        static_assert( (bindings::is_mutable< VectorTAU >::value) );
        assert( bindings::size(tau) >= bindings::size_row(a) );
        assert( bindings::size(work.select(real_type())) >=
                min_size_work( bindings::size_row(a) ));
        assert( bindings::size_column(a) >= 0 );
        assert( bindings::size_column(a) >= 0 );
        assert( bindings::size_minor(a) == 1 ||
                bindings::stride_minor(a) == 1 );
        assert( bindings::size_row(a) >= 0 );
        assert( bindings::stride_major(a) >= std::max< std::ptrdiff_t >(1,
                bindings::size_row(a)) );
        return detail::latrz( bindings::size_row(a), bindings::size_column(a),
                bindings::size_column(a), bindings::begin_value(a),
                bindings::stride_major(a), bindings::begin_value(tau),
                bindings::begin_value(work.select(real_type())) );
    }

    //
    // Static member function that
    // * Figures out the minimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member function
    // * Enables the unblocked algorithm (BLAS level 2)
    //
    template< typename MatrixA, typename VectorTAU >
    static std::ptrdiff_t invoke( MatrixA& a, VectorTAU& tau,
            minimal_workspace work ) {
        namespace bindings = ::boost::numeric::bindings;
        bindings::detail::array< real_type > tmp_work( min_size_work(
                bindings::size_row(a) ) );
        return invoke( a, tau, workspace( tmp_work ) );
    }

    //
    // Static member function that
    // * Figures out the optimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member
    // * Enables the blocked algorithm (BLAS level 3)
    //
    template< typename MatrixA, typename VectorTAU >
    static std::ptrdiff_t invoke( MatrixA& a, VectorTAU& tau,
            optimal_workspace work ) {
        namespace bindings = ::boost::numeric::bindings;
        return invoke( a, tau, minimal_workspace() );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array work.
    //
    static std::ptrdiff_t min_size_work( const std::ptrdiff_t m ) {
        return m;
    }
};

//
// This implementation is enabled if Value is a complex type.
//
template< typename Value >
struct latrz_impl< Value, typename std::enable_if< is_complex< Value >::value >::type > {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function for user-defined workspaces, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename MatrixA, typename VectorTAU, typename WORK >
    static std::ptrdiff_t invoke( MatrixA& a, VectorTAU& tau,
            detail::workspace1< WORK > work ) {
        namespace bindings = ::boost::numeric::bindings;
        static_assert( (bindings::is_column_major< MatrixA >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorTAU >::type >::type >::value) );
        static_assert( (bindings::is_mutable< MatrixA >::value) );
        static_assert( (bindings::is_mutable< VectorTAU >::value) );
        assert( bindings::size(tau) >= bindings::size_row(a) );
        assert( bindings::size(work.select(value_type())) >=
                min_size_work( bindings::size_row(a) ));
        assert( bindings::size_column(a) >= 0 );
        assert( bindings::size_column(a) >= 0 );
        assert( bindings::size_minor(a) == 1 ||
                bindings::stride_minor(a) == 1 );
        assert( bindings::size_row(a) >= 0 );
        assert( bindings::stride_major(a) >= std::max< std::ptrdiff_t >(1,
                bindings::size_row(a)) );
        return detail::latrz( bindings::size_row(a), bindings::size_column(a),
                bindings::size_column(a), bindings::begin_value(a),
                bindings::stride_major(a), bindings::begin_value(tau),
                bindings::begin_value(work.select(value_type())) );
    }

    //
    // Static member function that
    // * Figures out the minimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member function
    // * Enables the unblocked algorithm (BLAS level 2)
    //
    template< typename MatrixA, typename VectorTAU >
    static std::ptrdiff_t invoke( MatrixA& a, VectorTAU& tau,
            minimal_workspace work ) {
        namespace bindings = ::boost::numeric::bindings;
        bindings::detail::array< value_type > tmp_work( min_size_work(
                bindings::size_row(a) ) );
        return invoke( a, tau, workspace( tmp_work ) );
    }

    //
    // Static member function that
    // * Figures out the optimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member
    // * Enables the blocked algorithm (BLAS level 3)
    //
    template< typename MatrixA, typename VectorTAU >
    static std::ptrdiff_t invoke( MatrixA& a, VectorTAU& tau,
            optimal_workspace work ) {
        namespace bindings = ::boost::numeric::bindings;
        return invoke( a, tau, minimal_workspace() );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array work.
    //
    static std::ptrdiff_t min_size_work( const std::ptrdiff_t m ) {
        return m;
    }
};


//
// Functions for direct use. These functions are overloaded for temporaries,
// so that wrapped types can still be passed and used for write-access. In
// addition, if applicable, they are overloaded for user-defined workspaces.
// Calls to these functions are passed to the latrz_impl classes. In the 
// documentation, most overloads are collapsed to avoid a large number of
// prototypes which are very similar.
//

//
// Overloaded function for latrz. Its overload differs for
// * MatrixA&
// * VectorTAU&
// * User-defined workspace
//
template< typename MatrixA, typename VectorTAU, typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
latrz( MatrixA& a, VectorTAU& tau, Workspace work ) {
    return latrz_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( a, tau, work );
}

//
// Overloaded function for latrz. Its overload differs for
// * MatrixA&
// * VectorTAU&
// * Default workspace-type (optimal)
//
template< typename MatrixA, typename VectorTAU >
inline typename std::enable_if<! detail::is_workspace< VectorTAU >::value,
        std::ptrdiff_t >::type
latrz( MatrixA& a, VectorTAU& tau ) {
    return latrz_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( a, tau, optimal_workspace() );
}

//
// Overloaded function for latrz. Its overload differs for
// * const MatrixA&
// * VectorTAU&
// * User-defined workspace
//
template< typename MatrixA, typename VectorTAU, typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
latrz( const MatrixA& a, VectorTAU& tau, Workspace work ) {
    return latrz_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( a, tau, work );
}

//
// Overloaded function for latrz. Its overload differs for
// * const MatrixA&
// * VectorTAU&
// * Default workspace-type (optimal)
//
template< typename MatrixA, typename VectorTAU >
inline typename std::enable_if<! detail::is_workspace< VectorTAU >::value,
        std::ptrdiff_t >::type
latrz( const MatrixA& a, VectorTAU& tau ) {
    return latrz_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( a, tau, optimal_workspace() );
}

//
// Overloaded function for latrz. Its overload differs for
// * MatrixA&
// * const VectorTAU&
// * User-defined workspace
//
template< typename MatrixA, typename VectorTAU, typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
latrz( MatrixA& a, const VectorTAU& tau, Workspace work ) {
    return latrz_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( a, tau, work );
}

//
// Overloaded function for latrz. Its overload differs for
// * MatrixA&
// * const VectorTAU&
// * Default workspace-type (optimal)
//
template< typename MatrixA, typename VectorTAU >
inline typename std::enable_if<! detail::is_workspace< VectorTAU >::value,
        std::ptrdiff_t >::type
latrz( MatrixA& a, const VectorTAU& tau ) {
    return latrz_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( a, tau, optimal_workspace() );
}

//
// Overloaded function for latrz. Its overload differs for
// * const MatrixA&
// * const VectorTAU&
// * User-defined workspace
//
template< typename MatrixA, typename VectorTAU, typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
latrz( const MatrixA& a, const VectorTAU& tau, Workspace work ) {
    return latrz_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( a, tau, work );
}

//
// Overloaded function for latrz. Its overload differs for
// * const MatrixA&
// * const VectorTAU&
// * Default workspace-type (optimal)
//
template< typename MatrixA, typename VectorTAU >
inline typename std::enable_if<! detail::is_workspace< VectorTAU >::value,
        std::ptrdiff_t >::type
latrz( const MatrixA& a, const VectorTAU& tau ) {
    return latrz_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( a, tau, optimal_workspace() );
}

} // namespace lapack
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif