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

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_COMPUTATIONAL_STEIN_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_COMPUTATIONAL_STEIN_HPP

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
// The LAPACK-backend for stein is the netlib-compatible backend.
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
inline std::ptrdiff_t stein( const fortran_int_t n, const float* d,
        const float* e, const fortran_int_t m, const float* w,
        const fortran_int_t* iblock, const fortran_int_t* isplit, float* z,
        const fortran_int_t ldz, float* work, fortran_int_t* iwork,
        fortran_int_t* ifail ) {
    fortran_int_t info(0);
    LAPACK_SSTEIN( &n, d, e, &m, w, iblock, isplit, z, &ldz, work, iwork,
            ifail, &info );
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * double value-type.
//
inline std::ptrdiff_t stein( const fortran_int_t n, const double* d,
        const double* e, const fortran_int_t m, const double* w,
        const fortran_int_t* iblock, const fortran_int_t* isplit, double* z,
        const fortran_int_t ldz, double* work, fortran_int_t* iwork,
        fortran_int_t* ifail ) {
    fortran_int_t info(0);
    LAPACK_DSTEIN( &n, d, e, &m, w, iblock, isplit, z, &ldz, work, iwork,
            ifail, &info );
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<float> value-type.
//
inline std::ptrdiff_t stein( const fortran_int_t n, const float* d,
        const float* e, const fortran_int_t m, const float* w,
        const fortran_int_t* iblock, const fortran_int_t* isplit,
        std::complex<float>* z, const fortran_int_t ldz, float* work,
        fortran_int_t* iwork, fortran_int_t* ifail ) {
    fortran_int_t info(0);
    LAPACK_CSTEIN( &n, d, e, &m, w, iblock, isplit, z, &ldz, work, iwork,
            ifail, &info );
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<double> value-type.
//
inline std::ptrdiff_t stein( const fortran_int_t n, const double* d,
        const double* e, const fortran_int_t m, const double* w,
        const fortran_int_t* iblock, const fortran_int_t* isplit,
        std::complex<double>* z, const fortran_int_t ldz, double* work,
        fortran_int_t* iwork, fortran_int_t* ifail ) {
    fortran_int_t info(0);
    LAPACK_ZSTEIN( &n, d, e, &m, w, iblock, isplit, z, &ldz, work, iwork,
            ifail, &info );
    return info;
}

} // namespace detail

//
// Value-type based template class. Use this class if you need a type
// for dispatching to stein.
//
template< typename Value, typename Enable = void >
struct stein_impl {};

//
// This implementation is enabled if Value is a real type.
//
template< typename Value >
struct stein_impl< Value, typename std::enable_if< is_real< Value >::value >::type > {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function for user-defined workspaces, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename VectorD, typename VectorE, typename VectorW,
            typename VectorIBLOCK, typename VectorISPLIT, typename MatrixZ,
            typename VectorIFAIL, typename WORK, typename IWORK >
    static std::ptrdiff_t invoke( const fortran_int_t n, const VectorD& d,
            const VectorE& e, const fortran_int_t m, const VectorW& w,
            const VectorIBLOCK& iblock, const VectorISPLIT& isplit,
            MatrixZ& z, VectorIFAIL& ifail, detail::workspace2< WORK,
            IWORK > work ) {
        namespace bindings = ::boost::numeric::bindings;
        static_assert( (bindings::is_column_major< MatrixZ >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< VectorD >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorE >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< VectorD >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorW >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< VectorD >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixZ >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< VectorIBLOCK >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorISPLIT >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< VectorIBLOCK >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorIFAIL >::type >::type >::value) );
        static_assert( (bindings::is_mutable< MatrixZ >::value) );
        static_assert( (bindings::is_mutable< VectorIFAIL >::value) );
        assert( bindings::size(d) >= n );
        assert( bindings::size(e) >= n-1 );
        assert( bindings::size(isplit) >= n );
        assert( bindings::size(w) >= n );
        assert( bindings::size(work.select(fortran_int_t())) >=
                min_size_iwork( n ));
        assert( bindings::size(work.select(real_type())) >=
                min_size_work( n ));
        assert( bindings::size_minor(z) == 1 ||
                bindings::stride_minor(z) == 1 );
        assert( bindings::stride_major(z) >= std::max< std::ptrdiff_t >(1,
                n) );
        assert( n >= 0 );
        return detail::stein( n, bindings::begin_value(d),
                bindings::begin_value(e), m, bindings::begin_value(w),
                bindings::begin_value(iblock), bindings::begin_value(isplit),
                bindings::begin_value(z), bindings::stride_major(z),
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
    template< typename VectorD, typename VectorE, typename VectorW,
            typename VectorIBLOCK, typename VectorISPLIT, typename MatrixZ,
            typename VectorIFAIL >
    static std::ptrdiff_t invoke( const fortran_int_t n, const VectorD& d,
            const VectorE& e, const fortran_int_t m, const VectorW& w,
            const VectorIBLOCK& iblock, const VectorISPLIT& isplit,
            MatrixZ& z, VectorIFAIL& ifail, minimal_workspace work ) {
        namespace bindings = ::boost::numeric::bindings;
        bindings::detail::array< real_type > tmp_work( min_size_work( n ) );
        bindings::detail::array< fortran_int_t > tmp_iwork(
                min_size_iwork( n ) );
        return invoke( n, d, e, m, w, iblock, isplit, z, ifail,
                workspace( tmp_work, tmp_iwork ) );
    }

    //
    // Static member function that
    // * Figures out the optimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member
    // * Enables the blocked algorithm (BLAS level 3)
    //
    template< typename VectorD, typename VectorE, typename VectorW,
            typename VectorIBLOCK, typename VectorISPLIT, typename MatrixZ,
            typename VectorIFAIL >
    static std::ptrdiff_t invoke( const fortran_int_t n, const VectorD& d,
            const VectorE& e, const fortran_int_t m, const VectorW& w,
            const VectorIBLOCK& iblock, const VectorISPLIT& isplit,
            MatrixZ& z, VectorIFAIL& ifail, optimal_workspace work ) {
        namespace bindings = ::boost::numeric::bindings;
        return invoke( n, d, e, m, w, iblock, isplit, z, ifail,
                minimal_workspace() );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array work.
    //
    static std::ptrdiff_t min_size_work( const std::ptrdiff_t n ) {
        return 5*n;
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
struct stein_impl< Value, typename std::enable_if< is_complex< Value >::value >::type > {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function for user-defined workspaces, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename VectorD, typename VectorE, typename VectorW,
            typename VectorIBLOCK, typename VectorISPLIT, typename MatrixZ,
            typename VectorIFAIL, typename WORK, typename IWORK >
    static std::ptrdiff_t invoke( const fortran_int_t n, const VectorD& d,
            const VectorE& e, const fortran_int_t m, const VectorW& w,
            const VectorIBLOCK& iblock, const VectorISPLIT& isplit,
            MatrixZ& z, VectorIFAIL& ifail, detail::workspace2< WORK,
            IWORK > work ) {
        namespace bindings = ::boost::numeric::bindings;
        static_assert( (bindings::is_column_major< MatrixZ >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< VectorD >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorE >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< VectorD >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorW >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< VectorIBLOCK >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorISPLIT >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< VectorIBLOCK >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorIFAIL >::type >::type >::value) );
        static_assert( (bindings::is_mutable< MatrixZ >::value) );
        static_assert( (bindings::is_mutable< VectorIFAIL >::value) );
        assert( bindings::size(d) >= n );
        assert( bindings::size(e) >= n-1 );
        assert( bindings::size(isplit) >= n );
        assert( bindings::size(w) >= n );
        assert( bindings::size(work.select(fortran_int_t())) >=
                min_size_iwork( n ));
        assert( bindings::size(work.select(real_type())) >=
                min_size_work( n ));
        assert( bindings::size_minor(z) == 1 ||
                bindings::stride_minor(z) == 1 );
        assert( bindings::stride_major(z) >= std::max< std::ptrdiff_t >(1,
                n) );
        assert( n >= 0 );
        return detail::stein( n, bindings::begin_value(d),
                bindings::begin_value(e), m, bindings::begin_value(w),
                bindings::begin_value(iblock), bindings::begin_value(isplit),
                bindings::begin_value(z), bindings::stride_major(z),
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
    template< typename VectorD, typename VectorE, typename VectorW,
            typename VectorIBLOCK, typename VectorISPLIT, typename MatrixZ,
            typename VectorIFAIL >
    static std::ptrdiff_t invoke( const fortran_int_t n, const VectorD& d,
            const VectorE& e, const fortran_int_t m, const VectorW& w,
            const VectorIBLOCK& iblock, const VectorISPLIT& isplit,
            MatrixZ& z, VectorIFAIL& ifail, minimal_workspace work ) {
        namespace bindings = ::boost::numeric::bindings;
        bindings::detail::array< real_type > tmp_work( min_size_work( n ) );
        bindings::detail::array< fortran_int_t > tmp_iwork(
                min_size_iwork( n ) );
        return invoke( n, d, e, m, w, iblock, isplit, z, ifail,
                workspace( tmp_work, tmp_iwork ) );
    }

    //
    // Static member function that
    // * Figures out the optimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member
    // * Enables the blocked algorithm (BLAS level 3)
    //
    template< typename VectorD, typename VectorE, typename VectorW,
            typename VectorIBLOCK, typename VectorISPLIT, typename MatrixZ,
            typename VectorIFAIL >
    static std::ptrdiff_t invoke( const fortran_int_t n, const VectorD& d,
            const VectorE& e, const fortran_int_t m, const VectorW& w,
            const VectorIBLOCK& iblock, const VectorISPLIT& isplit,
            MatrixZ& z, VectorIFAIL& ifail, optimal_workspace work ) {
        namespace bindings = ::boost::numeric::bindings;
        return invoke( n, d, e, m, w, iblock, isplit, z, ifail,
                minimal_workspace() );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array work.
    //
    static std::ptrdiff_t min_size_work( const std::ptrdiff_t n ) {
        return 5*n;
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
// Functions for direct use. These functions are overloaded for temporaries,
// so that wrapped types can still be passed and used for write-access. In
// addition, if applicable, they are overloaded for user-defined workspaces.
// Calls to these functions are passed to the stein_impl classes. In the 
// documentation, most overloads are collapsed to avoid a large number of
// prototypes which are very similar.
//

//
// Overloaded function for stein. Its overload differs for
// * MatrixZ&
// * VectorIFAIL&
// * User-defined workspace
//
template< typename VectorD, typename VectorE, typename VectorW,
        typename VectorIBLOCK, typename VectorISPLIT, typename MatrixZ,
        typename VectorIFAIL, typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
stein( const fortran_int_t n, const VectorD& d, const VectorE& e,
        const fortran_int_t m, const VectorW& w,
        const VectorIBLOCK& iblock, const VectorISPLIT& isplit, MatrixZ& z,
        VectorIFAIL& ifail, Workspace work ) {
    return stein_impl< typename bindings::value_type<
            MatrixZ >::type >::invoke( n, d, e, m, w, iblock, isplit, z,
            ifail, work );
}

//
// Overloaded function for stein. Its overload differs for
// * MatrixZ&
// * VectorIFAIL&
// * Default workspace-type (optimal)
//
template< typename VectorD, typename VectorE, typename VectorW,
        typename VectorIBLOCK, typename VectorISPLIT, typename MatrixZ,
        typename VectorIFAIL >
inline typename std::enable_if<! detail::is_workspace< VectorIFAIL >::value,
        std::ptrdiff_t >::type
stein( const fortran_int_t n, const VectorD& d, const VectorE& e,
        const fortran_int_t m, const VectorW& w,
        const VectorIBLOCK& iblock, const VectorISPLIT& isplit, MatrixZ& z,
        VectorIFAIL& ifail ) {
    return stein_impl< typename bindings::value_type<
            MatrixZ >::type >::invoke( n, d, e, m, w, iblock, isplit, z,
            ifail, optimal_workspace() );
}

//
// Overloaded function for stein. Its overload differs for
// * const MatrixZ&
// * VectorIFAIL&
// * User-defined workspace
//
template< typename VectorD, typename VectorE, typename VectorW,
        typename VectorIBLOCK, typename VectorISPLIT, typename MatrixZ,
        typename VectorIFAIL, typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
stein( const fortran_int_t n, const VectorD& d, const VectorE& e,
        const fortran_int_t m, const VectorW& w,
        const VectorIBLOCK& iblock, const VectorISPLIT& isplit,
        const MatrixZ& z, VectorIFAIL& ifail, Workspace work ) {
    return stein_impl< typename bindings::value_type<
            MatrixZ >::type >::invoke( n, d, e, m, w, iblock, isplit, z,
            ifail, work );
}

//
// Overloaded function for stein. Its overload differs for
// * const MatrixZ&
// * VectorIFAIL&
// * Default workspace-type (optimal)
//
template< typename VectorD, typename VectorE, typename VectorW,
        typename VectorIBLOCK, typename VectorISPLIT, typename MatrixZ,
        typename VectorIFAIL >
inline typename std::enable_if<! detail::is_workspace< VectorIFAIL >::value,
        std::ptrdiff_t >::type
stein( const fortran_int_t n, const VectorD& d, const VectorE& e,
        const fortran_int_t m, const VectorW& w,
        const VectorIBLOCK& iblock, const VectorISPLIT& isplit,
        const MatrixZ& z, VectorIFAIL& ifail ) {
    return stein_impl< typename bindings::value_type<
            MatrixZ >::type >::invoke( n, d, e, m, w, iblock, isplit, z,
            ifail, optimal_workspace() );
}

//
// Overloaded function for stein. Its overload differs for
// * MatrixZ&
// * const VectorIFAIL&
// * User-defined workspace
//
template< typename VectorD, typename VectorE, typename VectorW,
        typename VectorIBLOCK, typename VectorISPLIT, typename MatrixZ,
        typename VectorIFAIL, typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
stein( const fortran_int_t n, const VectorD& d, const VectorE& e,
        const fortran_int_t m, const VectorW& w,
        const VectorIBLOCK& iblock, const VectorISPLIT& isplit, MatrixZ& z,
        const VectorIFAIL& ifail, Workspace work ) {
    return stein_impl< typename bindings::value_type<
            MatrixZ >::type >::invoke( n, d, e, m, w, iblock, isplit, z,
            ifail, work );
}

//
// Overloaded function for stein. Its overload differs for
// * MatrixZ&
// * const VectorIFAIL&
// * Default workspace-type (optimal)
//
template< typename VectorD, typename VectorE, typename VectorW,
        typename VectorIBLOCK, typename VectorISPLIT, typename MatrixZ,
        typename VectorIFAIL >
inline typename std::enable_if<! detail::is_workspace< VectorIFAIL >::value,
        std::ptrdiff_t >::type
stein( const fortran_int_t n, const VectorD& d, const VectorE& e,
        const fortran_int_t m, const VectorW& w,
        const VectorIBLOCK& iblock, const VectorISPLIT& isplit, MatrixZ& z,
        const VectorIFAIL& ifail ) {
    return stein_impl< typename bindings::value_type<
            MatrixZ >::type >::invoke( n, d, e, m, w, iblock, isplit, z,
            ifail, optimal_workspace() );
}

//
// Overloaded function for stein. Its overload differs for
// * const MatrixZ&
// * const VectorIFAIL&
// * User-defined workspace
//
template< typename VectorD, typename VectorE, typename VectorW,
        typename VectorIBLOCK, typename VectorISPLIT, typename MatrixZ,
        typename VectorIFAIL, typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
stein( const fortran_int_t n, const VectorD& d, const VectorE& e,
        const fortran_int_t m, const VectorW& w,
        const VectorIBLOCK& iblock, const VectorISPLIT& isplit,
        const MatrixZ& z, const VectorIFAIL& ifail, Workspace work ) {
    return stein_impl< typename bindings::value_type<
            MatrixZ >::type >::invoke( n, d, e, m, w, iblock, isplit, z,
            ifail, work );
}

//
// Overloaded function for stein. Its overload differs for
// * const MatrixZ&
// * const VectorIFAIL&
// * Default workspace-type (optimal)
//
template< typename VectorD, typename VectorE, typename VectorW,
        typename VectorIBLOCK, typename VectorISPLIT, typename MatrixZ,
        typename VectorIFAIL >
inline typename std::enable_if<! detail::is_workspace< VectorIFAIL >::value,
        std::ptrdiff_t >::type
stein( const fortran_int_t n, const VectorD& d, const VectorE& e,
        const fortran_int_t m, const VectorW& w,
        const VectorIBLOCK& iblock, const VectorISPLIT& isplit,
        const MatrixZ& z, const VectorIFAIL& ifail ) {
    return stein_impl< typename bindings::value_type<
            MatrixZ >::type >::invoke( n, d, e, m, w, iblock, isplit, z,
            ifail, optimal_workspace() );
}

} // namespace lapack
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
