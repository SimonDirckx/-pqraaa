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

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_AUXILIARY_LARFB_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_AUXILIARY_LARFB_HPP

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
// The LAPACK-backend for larfb is the netlib-compatible backend.
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
template< typename Side, typename Trans >
inline std::ptrdiff_t larfb( const Side side, const Trans trans,
        const char direct, const char storev, const fortran_int_t m,
        const fortran_int_t n, const fortran_int_t k, const float* v,
        const fortran_int_t ldv, const float* t, const fortran_int_t ldt,
        float* c, const fortran_int_t ldc, float* work,
        const fortran_int_t ldwork ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_SLARFB( &lapack_option< Side >::value, &lapack_option<
            Trans >::value, &direct, &storev, &m, &n, &k, v, &ldv, t, &ldt, c,
            &ldc, work, &ldwork );
#else
    LAPACK_SLARFB( &lapack_option< Side >::value, &lapack_option<
            Trans >::value, &direct, &storev, &m, &n, &k, v, &ldv, t, &ldt, c,
            &ldc, work, &ldwork ,1 ,1 ,1 ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * double value-type.
//
template< typename Side, typename Trans >
inline std::ptrdiff_t larfb( const Side side, const Trans trans,
        const char direct, const char storev, const fortran_int_t m,
        const fortran_int_t n, const fortran_int_t k, const double* v,
        const fortran_int_t ldv, const double* t, const fortran_int_t ldt,
        double* c, const fortran_int_t ldc, double* work,
        const fortran_int_t ldwork ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_DLARFB( &lapack_option< Side >::value, &lapack_option<
            Trans >::value, &direct, &storev, &m, &n, &k, v, &ldv, t, &ldt, c,
            &ldc, work, &ldwork );
#else
    LAPACK_DLARFB( &lapack_option< Side >::value, &lapack_option<
            Trans >::value, &direct, &storev, &m, &n, &k, v, &ldv, t, &ldt, c,
            &ldc, work, &ldwork ,1 ,1 ,1 ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<float> value-type.
//
template< typename Side, typename Trans >
inline std::ptrdiff_t larfb( const Side side, const Trans trans,
        const char direct, const char storev, const fortran_int_t m,
        const fortran_int_t n, const fortran_int_t k,
        const std::complex<float>* v, const fortran_int_t ldv,
        const std::complex<float>* t, const fortran_int_t ldt,
        std::complex<float>* c, const fortran_int_t ldc,
        std::complex<float>* work, const fortran_int_t ldwork ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_CLARFB( &lapack_option< Side >::value, &lapack_option<
            Trans >::value, &direct, &storev, &m, &n, &k, v, &ldv, t, &ldt, c,
            &ldc, work, &ldwork );
#else
    LAPACK_CLARFB( &lapack_option< Side >::value, &lapack_option<
            Trans >::value, &direct, &storev, &m, &n, &k, v, &ldv, t, &ldt, c,
            &ldc, work, &ldwork ,1 ,1 ,1 ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<double> value-type.
//
template< typename Side, typename Trans >
inline std::ptrdiff_t larfb( const Side side, const Trans trans,
        const char direct, const char storev, const fortran_int_t m,
        const fortran_int_t n, const fortran_int_t k,
        const std::complex<double>* v, const fortran_int_t ldv,
        const std::complex<double>* t, const fortran_int_t ldt,
        std::complex<double>* c, const fortran_int_t ldc,
        std::complex<double>* work, const fortran_int_t ldwork ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_ZLARFB( &lapack_option< Side >::value, &lapack_option<
            Trans >::value, &direct, &storev, &m, &n, &k, v, &ldv, t, &ldt, c,
            &ldc, work, &ldwork );
#else
    LAPACK_ZLARFB( &lapack_option< Side >::value, &lapack_option<
            Trans >::value, &direct, &storev, &m, &n, &k, v, &ldv, t, &ldt, c,
            &ldc, work, &ldwork ,1 ,1 ,1 ,1 );
#endif
    return info;
}

} // namespace detail

//
// Value-type based template class. Use this class if you need a type
// for dispatching to larfb.
//
template< typename Value, typename Enable = void >
struct larfb_impl {};

//
// This implementation is enabled if Value is a real type.
//
template< typename Value >
struct larfb_impl< Value, typename std::enable_if< is_real< Value >::value >::type > {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function for user-defined workspaces, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename Side, typename MatrixV, typename MatrixT,
            typename MatrixC, typename WORK >
    static std::ptrdiff_t invoke( const Side side, const char direct,
            const char storev, const MatrixV& v, const MatrixT& t, MatrixC& c,
            detail::workspace1< WORK > work ) {
        namespace bindings = ::boost::numeric::bindings;
        static_assert( (bindings::is_column_major< MatrixV >::value) );
        static_assert( (bindings::is_column_major< MatrixT >::value) );
        static_assert( (bindings::is_column_major< MatrixC >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixV >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixT >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixV >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixC >::type >::type >::value) );
        static_assert( (bindings::is_mutable< MatrixC >::value) );
        assert( bindings::size(work.select(real_type())) >=
                min_size_work( $CALL_MIN_SIZE ));
        assert( bindings::size_minor(c) == 1 ||
                bindings::stride_minor(c) == 1 );
        assert( bindings::size_minor(t) == 1 ||
                bindings::stride_minor(t) == 1 );
        assert( bindings::size_minor(v) == 1 ||
                bindings::stride_minor(v) == 1 );
        assert( bindings::stride_major(c) >= std::max< std::ptrdiff_t >(1,
                bindings::size_row(c)) );
        assert( bindings::stride_major(t) >= bindings::size_column(t) );
        assert( direct == 'F' || direct == 'B' );
        assert( storev == 'C' || storev == 'R' );
        return detail::larfb( side, trans(), direct, storev,
                bindings::size_row(c), bindings::size_column(c),
                bindings::size_column(t), bindings::begin_value(v),
                bindings::stride_major(v), bindings::begin_value(t),
                bindings::stride_major(t), bindings::begin_value(c),
                bindings::stride_major(c), bindings::begin_value(work),
                bindings::stride_major(work) );
    }

    //
    // Static member function that
    // * Figures out the minimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member function
    // * Enables the unblocked algorithm (BLAS level 2)
    //
    template< typename Side, typename MatrixV, typename MatrixT,
            typename MatrixC >
    static std::ptrdiff_t invoke( const Side side, const char direct,
            const char storev, const MatrixV& v, const MatrixT& t, MatrixC& c,
            minimal_workspace work ) {
        namespace bindings = ::boost::numeric::bindings;
        bindings::detail::array< real_type > tmp_work( min_size_work(
                $CALL_MIN_SIZE ) );
        return invoke( side, direct, storev, v, t, c, workspace( tmp_work ) );
    }

    //
    // Static member function that
    // * Figures out the optimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member
    // * Enables the blocked algorithm (BLAS level 3)
    //
    template< typename Side, typename MatrixV, typename MatrixT,
            typename MatrixC >
    static std::ptrdiff_t invoke( const Side side, const char direct,
            const char storev, const MatrixV& v, const MatrixT& t, MatrixC& c,
            optimal_workspace work ) {
        namespace bindings = ::boost::numeric::bindings;
        return invoke( side, direct, storev, v, t, c, minimal_workspace() );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array work.
    //
    template< $TYPES >
    static std::ptrdiff_t min_size_work( $ARGUMENTS ) {
        $MIN_SIZE_IMPLEMENTATION
    }
};

//
// This implementation is enabled if Value is a complex type.
//
template< typename Value >
struct larfb_impl< Value, typename std::enable_if< is_complex< Value >::value >::type > {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function for user-defined workspaces, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename Side, typename MatrixV, typename MatrixT,
            typename MatrixC, typename WORK >
    static std::ptrdiff_t invoke( const Side side, const char direct,
            const char storev, const MatrixV& v, const MatrixT& t, MatrixC& c,
            detail::workspace1< WORK > work ) {
        namespace bindings = ::boost::numeric::bindings;
        static_assert( (bindings::is_column_major< MatrixV >::value) );
        static_assert( (bindings::is_column_major< MatrixT >::value) );
        static_assert( (bindings::is_column_major< MatrixC >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixV >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixT >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixV >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixC >::type >::type >::value) );
        static_assert( (bindings::is_mutable< MatrixC >::value) );
        assert( bindings::size(work.select(value_type())) >=
                min_size_work( $CALL_MIN_SIZE ));
        assert( bindings::size_minor(c) == 1 ||
                bindings::stride_minor(c) == 1 );
        assert( bindings::size_minor(t) == 1 ||
                bindings::stride_minor(t) == 1 );
        assert( bindings::size_minor(v) == 1 ||
                bindings::stride_minor(v) == 1 );
        assert( bindings::stride_major(c) >= std::max< std::ptrdiff_t >(1,
                bindings::size_row(c)) );
        assert( bindings::stride_major(t) >= bindings::size_column(t) );
        assert( direct == 'F' || direct == 'B' );
        assert( storev == 'C' || storev == 'R' );
        return detail::larfb( side, trans(), direct, storev,
                bindings::size_row(c), bindings::size_column(c),
                bindings::size_column(t), bindings::begin_value(v),
                bindings::stride_major(v), bindings::begin_value(t),
                bindings::stride_major(t), bindings::begin_value(c),
                bindings::stride_major(c), bindings::begin_value(work),
                bindings::stride_major(work) );
    }

    //
    // Static member function that
    // * Figures out the minimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member function
    // * Enables the unblocked algorithm (BLAS level 2)
    //
    template< typename Side, typename MatrixV, typename MatrixT,
            typename MatrixC >
    static std::ptrdiff_t invoke( const Side side, const char direct,
            const char storev, const MatrixV& v, const MatrixT& t, MatrixC& c,
            minimal_workspace work ) {
        namespace bindings = ::boost::numeric::bindings;
        bindings::detail::array< value_type > tmp_work( min_size_work(
                $CALL_MIN_SIZE ) );
        return invoke( side, direct, storev, v, t, c, workspace( tmp_work ) );
    }

    //
    // Static member function that
    // * Figures out the optimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member
    // * Enables the blocked algorithm (BLAS level 3)
    //
    template< typename Side, typename MatrixV, typename MatrixT,
            typename MatrixC >
    static std::ptrdiff_t invoke( const Side side, const char direct,
            const char storev, const MatrixV& v, const MatrixT& t, MatrixC& c,
            optimal_workspace work ) {
        namespace bindings = ::boost::numeric::bindings;
        return invoke( side, direct, storev, v, t, c, minimal_workspace() );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array work.
    //
    template< $TYPES >
    static std::ptrdiff_t min_size_work( $ARGUMENTS ) {
        $MIN_SIZE_IMPLEMENTATION
    }
};


//
// Functions for direct use. These functions are overloaded for temporaries,
// so that wrapped types can still be passed and used for write-access. In
// addition, if applicable, they are overloaded for user-defined workspaces.
// Calls to these functions are passed to the larfb_impl classes. In the 
// documentation, most overloads are collapsed to avoid a large number of
// prototypes which are very similar.
//

//
// Overloaded function for larfb. Its overload differs for
// * MatrixC&
// * User-defined workspace
//
template< typename Side, typename MatrixV, typename MatrixT, typename MatrixC,
        typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
larfb( const Side side, const char direct, const char storev,
        const MatrixV& v, const MatrixT& t, MatrixC& c, Workspace work ) {
    return larfb_impl< typename bindings::value_type<
            MatrixV >::type >::invoke( side, direct, storev, v, t, c, work );
}

//
// Overloaded function for larfb. Its overload differs for
// * MatrixC&
// * Default workspace-type (optimal)
//
template< typename Side, typename MatrixV, typename MatrixT, typename MatrixC >
inline typename std::enable_if<! detail::is_workspace< MatrixC >::value,
        std::ptrdiff_t >::type
larfb( const Side side, const char direct, const char storev,
        const MatrixV& v, const MatrixT& t, MatrixC& c ) {
    return larfb_impl< typename bindings::value_type<
            MatrixV >::type >::invoke( side, direct, storev, v, t, c,
            optimal_workspace() );
}

//
// Overloaded function for larfb. Its overload differs for
// * const MatrixC&
// * User-defined workspace
//
template< typename Side, typename MatrixV, typename MatrixT, typename MatrixC,
        typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
larfb( const Side side, const char direct, const char storev,
        const MatrixV& v, const MatrixT& t, const MatrixC& c,
        Workspace work ) {
    return larfb_impl< typename bindings::value_type<
            MatrixV >::type >::invoke( side, direct, storev, v, t, c, work );
}

//
// Overloaded function for larfb. Its overload differs for
// * const MatrixC&
// * Default workspace-type (optimal)
//
template< typename Side, typename MatrixV, typename MatrixT, typename MatrixC >
inline typename std::enable_if<! detail::is_workspace< MatrixC >::value,
        std::ptrdiff_t >::type
larfb( const Side side, const char direct, const char storev,
        const MatrixV& v, const MatrixT& t, const MatrixC& c ) {
    return larfb_impl< typename bindings::value_type<
            MatrixV >::type >::invoke( side, direct, storev, v, t, c,
            optimal_workspace() );
}

} // namespace lapack
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
