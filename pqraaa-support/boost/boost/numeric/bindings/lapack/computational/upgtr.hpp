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

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_COMPUTATIONAL_UPGTR_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_COMPUTATIONAL_UPGTR_HPP

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
// The LAPACK-backend for upgtr is the netlib-compatible backend.
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
// * complex<float> value-type.
//
inline std::ptrdiff_t upgtr( const char uplo, const fortran_int_t n,
        const std::complex<float>* ap, const std::complex<float>* tau,
        std::complex<float>* q, const fortran_int_t ldq,
        std::complex<float>* work ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_CUPGTR( &uplo, &n, ap, tau, q, &ldq, work, &info );
#else
    LAPACK_CUPGTR( &uplo, &n, ap, tau, q, &ldq, work, &info ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<double> value-type.
//
inline std::ptrdiff_t upgtr( const char uplo, const fortran_int_t n,
        const std::complex<double>* ap, const std::complex<double>* tau,
        std::complex<double>* q, const fortran_int_t ldq,
        std::complex<double>* work ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_ZUPGTR( &uplo, &n, ap, tau, q, &ldq, work, &info );
#else
    LAPACK_ZUPGTR( &uplo, &n, ap, tau, q, &ldq, work, &info ,1 );
#endif
    return info;
}

} // namespace detail

//
// Value-type based template class. Use this class if you need a type
// for dispatching to upgtr.
//
template< typename Value >
struct upgtr_impl {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function for user-defined workspaces, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename VectorAP, typename VectorTAU, typename MatrixQ,
            typename WORK >
    static std::ptrdiff_t invoke( const char uplo, const VectorAP& ap,
            const VectorTAU& tau, MatrixQ& q, detail::workspace1<
            WORK > work ) {
        namespace bindings = ::boost::numeric::bindings;
        static_assert( (bindings::is_column_major< MatrixQ >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< VectorAP >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorTAU >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< VectorAP >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixQ >::type >::type >::value) );
        static_assert( (bindings::is_mutable< MatrixQ >::value) );
        assert( bindings::size(ap) >=
                bindings::size_column(q)*(bindings::size_column(q)+1)/2 );
        assert( bindings::size(tau) >= bindings::size_column(q)-1 );
        assert( bindings::size(work.select(value_type())) >=
                min_size_work( bindings::size_column(q) ));
        assert( bindings::size_column(q) >= 0 );
        assert( bindings::size_minor(q) == 1 ||
                bindings::stride_minor(q) == 1 );
        assert( bindings::stride_major(q) >= std::max< std::ptrdiff_t >(1,
                bindings::size_column(q)) );
        return detail::upgtr( uplo, bindings::size_column(q),
                bindings::begin_value(ap), bindings::begin_value(tau),
                bindings::begin_value(q), bindings::stride_major(q),
                bindings::begin_value(work.select(value_type())) );
    }

    //
    // Static member function that
    // * Figures out the minimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member function
    // * Enables the unblocked algorithm (BLAS level 2)
    //
    template< typename VectorAP, typename VectorTAU, typename MatrixQ >
    static std::ptrdiff_t invoke( const char uplo, const VectorAP& ap,
            const VectorTAU& tau, MatrixQ& q, minimal_workspace work ) {
        namespace bindings = ::boost::numeric::bindings;
        bindings::detail::array< value_type > tmp_work( min_size_work(
                bindings::size_column(q) ) );
        return invoke( uplo, ap, tau, q, workspace( tmp_work ) );
    }

    //
    // Static member function that
    // * Figures out the optimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member
    // * Enables the blocked algorithm (BLAS level 3)
    //
    template< typename VectorAP, typename VectorTAU, typename MatrixQ >
    static std::ptrdiff_t invoke( const char uplo, const VectorAP& ap,
            const VectorTAU& tau, MatrixQ& q, optimal_workspace work ) {
        namespace bindings = ::boost::numeric::bindings;
        return invoke( uplo, ap, tau, q, minimal_workspace() );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array work.
    //
    static std::ptrdiff_t min_size_work( const std::ptrdiff_t n ) {
        return n-1;
    }
};


//
// Functions for direct use. These functions are overloaded for temporaries,
// so that wrapped types can still be passed and used for write-access. In
// addition, if applicable, they are overloaded for user-defined workspaces.
// Calls to these functions are passed to the upgtr_impl classes. In the 
// documentation, most overloads are collapsed to avoid a large number of
// prototypes which are very similar.
//

//
// Overloaded function for upgtr. Its overload differs for
// * MatrixQ&
// * User-defined workspace
//
template< typename VectorAP, typename VectorTAU, typename MatrixQ,
        typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
upgtr( const char uplo, const VectorAP& ap, const VectorTAU& tau,
        MatrixQ& q, Workspace work ) {
    return upgtr_impl< typename bindings::value_type<
            VectorAP >::type >::invoke( uplo, ap, tau, q, work );
}

//
// Overloaded function for upgtr. Its overload differs for
// * MatrixQ&
// * Default workspace-type (optimal)
//
template< typename VectorAP, typename VectorTAU, typename MatrixQ >
inline typename std::enable_if<! detail::is_workspace< MatrixQ >::value,
        std::ptrdiff_t >::type
upgtr( const char uplo, const VectorAP& ap, const VectorTAU& tau,
        MatrixQ& q ) {
    return upgtr_impl< typename bindings::value_type<
            VectorAP >::type >::invoke( uplo, ap, tau, q,
            optimal_workspace() );
}

//
// Overloaded function for upgtr. Its overload differs for
// * const MatrixQ&
// * User-defined workspace
//
template< typename VectorAP, typename VectorTAU, typename MatrixQ,
        typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
upgtr( const char uplo, const VectorAP& ap, const VectorTAU& tau,
        const MatrixQ& q, Workspace work ) {
    return upgtr_impl< typename bindings::value_type<
            VectorAP >::type >::invoke( uplo, ap, tau, q, work );
}

//
// Overloaded function for upgtr. Its overload differs for
// * const MatrixQ&
// * Default workspace-type (optimal)
//
template< typename VectorAP, typename VectorTAU, typename MatrixQ >
inline typename std::enable_if<! detail::is_workspace< MatrixQ >::value,
        std::ptrdiff_t >::type
upgtr( const char uplo, const VectorAP& ap, const VectorTAU& tau,
        const MatrixQ& q ) {
    return upgtr_impl< typename bindings::value_type<
            VectorAP >::type >::invoke( uplo, ap, tau, q,
            optimal_workspace() );
}

} // namespace lapack
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
