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

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_COMPUTATIONAL_TRSNA_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_COMPUTATIONAL_TRSNA_HPP

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
// The LAPACK-backend for trsna is the netlib-compatible backend.
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
inline std::ptrdiff_t trsna( const char job, const char howmny,
        const logical_t* select, const fortran_int_t n, const float* t,
        const fortran_int_t ldt, const float* vl, const fortran_int_t ldvl,
        const float* vr, const fortran_int_t ldvr, float* s, float* sep,
        const fortran_int_t mm, fortran_int_t& m, float* work,
        const fortran_int_t ldwork, fortran_int_t* iwork ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_STRSNA( &job, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr,
            s, sep, &mm, &m, work, &ldwork, iwork, &info );
#else
    LAPACK_STRSNA( &job, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr,
            s, sep, &mm, &m, work, &ldwork, iwork, &info ,1 ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * double value-type.
//
inline std::ptrdiff_t trsna( const char job, const char howmny,
        const logical_t* select, const fortran_int_t n, const double* t,
        const fortran_int_t ldt, const double* vl, const fortran_int_t ldvl,
        const double* vr, const fortran_int_t ldvr, double* s, double* sep,
        const fortran_int_t mm, fortran_int_t& m, double* work,
        const fortran_int_t ldwork, fortran_int_t* iwork ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_DTRSNA( &job, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr,
            s, sep, &mm, &m, work, &ldwork, iwork, &info );
#else
    LAPACK_DTRSNA( &job, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr,
            s, sep, &mm, &m, work, &ldwork, iwork, &info ,1 ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<float> value-type.
//
inline std::ptrdiff_t trsna( const char job, const char howmny,
        const logical_t* select, const fortran_int_t n,
        const std::complex<float>* t, const fortran_int_t ldt,
        const std::complex<float>* vl, const fortran_int_t ldvl,
        const std::complex<float>* vr, const fortran_int_t ldvr, float* s,
        float* sep, const fortran_int_t mm, fortran_int_t& m,
        std::complex<float>* work, const fortran_int_t ldwork, float* rwork ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_CTRSNA( &job, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr,
            s, sep, &mm, &m, work, &ldwork, rwork, &info );
#else
    LAPACK_CTRSNA( &job, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr,
            s, sep, &mm, &m, work, &ldwork, rwork, &info ,1 ,1 );
#endif
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<double> value-type.
//
inline std::ptrdiff_t trsna( const char job, const char howmny,
        const logical_t* select, const fortran_int_t n,
        const std::complex<double>* t, const fortran_int_t ldt,
        const std::complex<double>* vl, const fortran_int_t ldvl,
        const std::complex<double>* vr, const fortran_int_t ldvr, double* s,
        double* sep, const fortran_int_t mm, fortran_int_t& m,
        std::complex<double>* work, const fortran_int_t ldwork,
        double* rwork ) {
    fortran_int_t info(0);
#ifndef LAPACK_FORTRAN_STRLEN_END
    LAPACK_ZTRSNA( &job, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr,
            s, sep, &mm, &m, work, &ldwork, rwork, &info );
#else
    LAPACK_ZTRSNA( &job, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr,
            s, sep, &mm, &m, work, &ldwork, rwork, &info ,1 ,1 );
#endif
    return info;
}

} // namespace detail

//
// Value-type based template class. Use this class if you need a type
// for dispatching to trsna.
//
template< typename Value, typename Enable = void >
struct trsna_impl {};

//
// This implementation is enabled if Value is a real type.
//
template< typename Value >
struct trsna_impl< Value, typename std::enable_if< is_real< Value >::value >::type > {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function for user-defined workspaces, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename VectorSELECT, typename MatrixT, typename MatrixVL,
            typename MatrixVR, typename VectorS, typename VectorSEP,
            typename WORK, typename IWORK >
    static std::ptrdiff_t invoke( const char job, const char howmny,
            const VectorSELECT& select, const MatrixT& t, const MatrixVL& vl,
            const MatrixVR& vr, VectorS& s, VectorSEP& sep,
            const fortran_int_t mm, fortran_int_t& m,
            detail::workspace2< WORK, IWORK > work ) {
        namespace bindings = ::boost::numeric::bindings;
        static_assert( (bindings::is_column_major< MatrixT >::value) );
        static_assert( (bindings::is_column_major< MatrixVL >::value) );
        static_assert( (bindings::is_column_major< MatrixVR >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixT >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixVL >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixT >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixVR >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixT >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorS >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixT >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorSEP >::type >::type >::value) );
        static_assert( (bindings::is_mutable< VectorS >::value) );
        static_assert( (bindings::is_mutable< VectorSEP >::value) );
        assert( bindings::size(work.select(fortran_int_t())) >=
                min_size_iwork( $CALL_MIN_SIZE ));
        assert( bindings::size(work.select(real_type())) >=
                min_size_work( $CALL_MIN_SIZE ));
        assert( bindings::size_column(t) >= 0 );
        assert( bindings::size_minor(t) == 1 ||
                bindings::stride_minor(t) == 1 );
        assert( bindings::size_minor(vl) == 1 ||
                bindings::stride_minor(vl) == 1 );
        assert( bindings::size_minor(vr) == 1 ||
                bindings::stride_minor(vr) == 1 );
        assert( bindings::stride_major(t) >= std::max< std::ptrdiff_t >(1,
                bindings::size_column(t)) );
        assert( howmny == 'A' || howmny == 'S' );
        assert( job == 'E' || job == 'V' || job == 'B' );
        return detail::trsna( job, howmny, bindings::begin_value(select),
                bindings::size_column(t), bindings::begin_value(t),
                bindings::stride_major(t), bindings::begin_value(vl),
                bindings::stride_major(vl), bindings::begin_value(vr),
                bindings::stride_major(vr), bindings::begin_value(s),
                bindings::begin_value(sep), mm, m,
                bindings::begin_value(work), bindings::stride_major(work),
                bindings::begin_value(work.select(fortran_int_t())) );
    }

    //
    // Static member function that
    // * Figures out the minimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member function
    // * Enables the unblocked algorithm (BLAS level 2)
    //
    template< typename VectorSELECT, typename MatrixT, typename MatrixVL,
            typename MatrixVR, typename VectorS, typename VectorSEP >
    static std::ptrdiff_t invoke( const char job, const char howmny,
            const VectorSELECT& select, const MatrixT& t, const MatrixVL& vl,
            const MatrixVR& vr, VectorS& s, VectorSEP& sep,
            const fortran_int_t mm, fortran_int_t& m,
            minimal_workspace work ) {
        namespace bindings = ::boost::numeric::bindings;
        bindings::detail::array< real_type > tmp_work( min_size_work(
                $CALL_MIN_SIZE ) );
        bindings::detail::array< fortran_int_t > tmp_iwork(
                min_size_iwork( $CALL_MIN_SIZE ) );
        return invoke( job, howmny, select, t, vl, vr, s, sep, mm, m,
                workspace( tmp_work, tmp_iwork ) );
    }

    //
    // Static member function that
    // * Figures out the optimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member
    // * Enables the blocked algorithm (BLAS level 3)
    //
    template< typename VectorSELECT, typename MatrixT, typename MatrixVL,
            typename MatrixVR, typename VectorS, typename VectorSEP >
    static std::ptrdiff_t invoke( const char job, const char howmny,
            const VectorSELECT& select, const MatrixT& t, const MatrixVL& vl,
            const MatrixVR& vr, VectorS& s, VectorSEP& sep,
            const fortran_int_t mm, fortran_int_t& m,
            optimal_workspace work ) {
        namespace bindings = ::boost::numeric::bindings;
        return invoke( job, howmny, select, t, vl, vr, s, sep, mm, m,
                minimal_workspace() );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array work.
    //
    template< $TYPES >
    static std::ptrdiff_t min_size_work( $ARGUMENTS ) {
        $MIN_SIZE_IMPLEMENTATION
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array iwork.
    //
    template< $TYPES >
    static std::ptrdiff_t min_size_iwork( $ARGUMENTS ) {
        $MIN_SIZE_IMPLEMENTATION
    }
};

//
// This implementation is enabled if Value is a complex type.
//
template< typename Value >
struct trsna_impl< Value, typename std::enable_if< is_complex< Value >::value >::type > {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function for user-defined workspaces, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename VectorSELECT, typename MatrixT, typename MatrixVL,
            typename MatrixVR, typename VectorS, typename VectorSEP,
            typename WORK, typename RWORK >
    static std::ptrdiff_t invoke( const char job, const char howmny,
            const VectorSELECT& select, const MatrixT& t, const MatrixVL& vl,
            const MatrixVR& vr, VectorS& s, VectorSEP& sep,
            const fortran_int_t mm, fortran_int_t& m,
            detail::workspace2< WORK, RWORK > work ) {
        namespace bindings = ::boost::numeric::bindings;
        static_assert( (bindings::is_column_major< MatrixT >::value) );
        static_assert( (bindings::is_column_major< MatrixVL >::value) );
        static_assert( (bindings::is_column_major< MatrixVR >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< VectorS >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorSEP >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixT >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixVL >::type >::type >::value) );
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< MatrixT >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                MatrixVR >::type >::type >::value) );
        static_assert( (bindings::is_mutable< VectorS >::value) );
        static_assert( (bindings::is_mutable< VectorSEP >::value) );
        assert( bindings::size(work.select(real_type())) >=
                min_size_rwork( $CALL_MIN_SIZE ));
        assert( bindings::size(work.select(value_type())) >=
                min_size_work( $CALL_MIN_SIZE ));
        assert( bindings::size_column(t) >= 0 );
        assert( bindings::size_minor(t) == 1 ||
                bindings::stride_minor(t) == 1 );
        assert( bindings::size_minor(vl) == 1 ||
                bindings::stride_minor(vl) == 1 );
        assert( bindings::size_minor(vr) == 1 ||
                bindings::stride_minor(vr) == 1 );
        assert( bindings::stride_major(t) >= std::max< std::ptrdiff_t >(1,
                bindings::size_column(t)) );
        assert( howmny == 'A' || howmny == 'S' );
        assert( job == 'E' || job == 'V' || job == 'B' );
        return detail::trsna( job, howmny, bindings::begin_value(select),
                bindings::size_column(t), bindings::begin_value(t),
                bindings::stride_major(t), bindings::begin_value(vl),
                bindings::stride_major(vl), bindings::begin_value(vr),
                bindings::stride_major(vr), bindings::begin_value(s),
                bindings::begin_value(sep), mm, m,
                bindings::begin_value(work), bindings::stride_major(work),
                bindings::begin_value(work.select(real_type())) );
    }

    //
    // Static member function that
    // * Figures out the minimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member function
    // * Enables the unblocked algorithm (BLAS level 2)
    //
    template< typename VectorSELECT, typename MatrixT, typename MatrixVL,
            typename MatrixVR, typename VectorS, typename VectorSEP >
    static std::ptrdiff_t invoke( const char job, const char howmny,
            const VectorSELECT& select, const MatrixT& t, const MatrixVL& vl,
            const MatrixVR& vr, VectorS& s, VectorSEP& sep,
            const fortran_int_t mm, fortran_int_t& m,
            minimal_workspace work ) {
        namespace bindings = ::boost::numeric::bindings;
        bindings::detail::array< value_type > tmp_work( min_size_work(
                $CALL_MIN_SIZE ) );
        bindings::detail::array< real_type > tmp_rwork( min_size_rwork(
                $CALL_MIN_SIZE ) );
        return invoke( job, howmny, select, t, vl, vr, s, sep, mm, m,
                workspace( tmp_work, tmp_rwork ) );
    }

    //
    // Static member function that
    // * Figures out the optimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member
    // * Enables the blocked algorithm (BLAS level 3)
    //
    template< typename VectorSELECT, typename MatrixT, typename MatrixVL,
            typename MatrixVR, typename VectorS, typename VectorSEP >
    static std::ptrdiff_t invoke( const char job, const char howmny,
            const VectorSELECT& select, const MatrixT& t, const MatrixVL& vl,
            const MatrixVR& vr, VectorS& s, VectorSEP& sep,
            const fortran_int_t mm, fortran_int_t& m,
            optimal_workspace work ) {
        namespace bindings = ::boost::numeric::bindings;
        return invoke( job, howmny, select, t, vl, vr, s, sep, mm, m,
                minimal_workspace() );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array work.
    //
    template< $TYPES >
    static std::ptrdiff_t min_size_work( $ARGUMENTS ) {
        $MIN_SIZE_IMPLEMENTATION
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array rwork.
    //
    template< $TYPES >
    static std::ptrdiff_t min_size_rwork( $ARGUMENTS ) {
        $MIN_SIZE_IMPLEMENTATION
    }
};


//
// Functions for direct use. These functions are overloaded for temporaries,
// so that wrapped types can still be passed and used for write-access. In
// addition, if applicable, they are overloaded for user-defined workspaces.
// Calls to these functions are passed to the trsna_impl classes. In the 
// documentation, most overloads are collapsed to avoid a large number of
// prototypes which are very similar.
//

//
// Overloaded function for trsna. Its overload differs for
// * VectorS&
// * VectorSEP&
// * User-defined workspace
//
template< typename VectorSELECT, typename MatrixT, typename MatrixVL,
        typename MatrixVR, typename VectorS, typename VectorSEP,
        typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
trsna( const char job, const char howmny, const VectorSELECT& select,
        const MatrixT& t, const MatrixVL& vl, const MatrixVR& vr, VectorS& s,
        VectorSEP& sep, const fortran_int_t mm, fortran_int_t& m,
        Workspace work ) {
    return trsna_impl< typename bindings::value_type<
            MatrixT >::type >::invoke( job, howmny, select, t, vl, vr, s, sep,
            mm, m, work );
}

//
// Overloaded function for trsna. Its overload differs for
// * VectorS&
// * VectorSEP&
// * Default workspace-type (optimal)
//
template< typename VectorSELECT, typename MatrixT, typename MatrixVL,
        typename MatrixVR, typename VectorS, typename VectorSEP >
inline typename std::enable_if<! detail::is_workspace< VectorSEP >::value,
        std::ptrdiff_t >::type
trsna( const char job, const char howmny, const VectorSELECT& select,
        const MatrixT& t, const MatrixVL& vl, const MatrixVR& vr, VectorS& s,
        VectorSEP& sep, const fortran_int_t mm, fortran_int_t& m ) {
    return trsna_impl< typename bindings::value_type<
            MatrixT >::type >::invoke( job, howmny, select, t, vl, vr, s, sep,
            mm, m, optimal_workspace() );
}

//
// Overloaded function for trsna. Its overload differs for
// * const VectorS&
// * VectorSEP&
// * User-defined workspace
//
template< typename VectorSELECT, typename MatrixT, typename MatrixVL,
        typename MatrixVR, typename VectorS, typename VectorSEP,
        typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
trsna( const char job, const char howmny, const VectorSELECT& select,
        const MatrixT& t, const MatrixVL& vl, const MatrixVR& vr,
        const VectorS& s, VectorSEP& sep, const fortran_int_t mm,
        fortran_int_t& m, Workspace work ) {
    return trsna_impl< typename bindings::value_type<
            MatrixT >::type >::invoke( job, howmny, select, t, vl, vr, s, sep,
            mm, m, work );
}

//
// Overloaded function for trsna. Its overload differs for
// * const VectorS&
// * VectorSEP&
// * Default workspace-type (optimal)
//
template< typename VectorSELECT, typename MatrixT, typename MatrixVL,
        typename MatrixVR, typename VectorS, typename VectorSEP >
inline typename std::enable_if<! detail::is_workspace< VectorSEP >::value,
        std::ptrdiff_t >::type
trsna( const char job, const char howmny, const VectorSELECT& select,
        const MatrixT& t, const MatrixVL& vl, const MatrixVR& vr,
        const VectorS& s, VectorSEP& sep, const fortran_int_t mm,
        fortran_int_t& m ) {
    return trsna_impl< typename bindings::value_type<
            MatrixT >::type >::invoke( job, howmny, select, t, vl, vr, s, sep,
            mm, m, optimal_workspace() );
}

//
// Overloaded function for trsna. Its overload differs for
// * VectorS&
// * const VectorSEP&
// * User-defined workspace
//
template< typename VectorSELECT, typename MatrixT, typename MatrixVL,
        typename MatrixVR, typename VectorS, typename VectorSEP,
        typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
trsna( const char job, const char howmny, const VectorSELECT& select,
        const MatrixT& t, const MatrixVL& vl, const MatrixVR& vr, VectorS& s,
        const VectorSEP& sep, const fortran_int_t mm,
        fortran_int_t& m, Workspace work ) {
    return trsna_impl< typename bindings::value_type<
            MatrixT >::type >::invoke( job, howmny, select, t, vl, vr, s, sep,
            mm, m, work );
}

//
// Overloaded function for trsna. Its overload differs for
// * VectorS&
// * const VectorSEP&
// * Default workspace-type (optimal)
//
template< typename VectorSELECT, typename MatrixT, typename MatrixVL,
        typename MatrixVR, typename VectorS, typename VectorSEP >
inline typename std::enable_if<! detail::is_workspace< VectorSEP >::value,
        std::ptrdiff_t >::type
trsna( const char job, const char howmny, const VectorSELECT& select,
        const MatrixT& t, const MatrixVL& vl, const MatrixVR& vr, VectorS& s,
        const VectorSEP& sep, const fortran_int_t mm,
        fortran_int_t& m ) {
    return trsna_impl< typename bindings::value_type<
            MatrixT >::type >::invoke( job, howmny, select, t, vl, vr, s, sep,
            mm, m, optimal_workspace() );
}

//
// Overloaded function for trsna. Its overload differs for
// * const VectorS&
// * const VectorSEP&
// * User-defined workspace
//
template< typename VectorSELECT, typename MatrixT, typename MatrixVL,
        typename MatrixVR, typename VectorS, typename VectorSEP,
        typename Workspace >
inline typename std::enable_if< detail::is_workspace< Workspace >::value,
        std::ptrdiff_t >::type
trsna( const char job, const char howmny, const VectorSELECT& select,
        const MatrixT& t, const MatrixVL& vl, const MatrixVR& vr,
        const VectorS& s, const VectorSEP& sep, const fortran_int_t mm,
        fortran_int_t& m, Workspace work ) {
    return trsna_impl< typename bindings::value_type<
            MatrixT >::type >::invoke( job, howmny, select, t, vl, vr, s, sep,
            mm, m, work );
}

//
// Overloaded function for trsna. Its overload differs for
// * const VectorS&
// * const VectorSEP&
// * Default workspace-type (optimal)
//
template< typename VectorSELECT, typename MatrixT, typename MatrixVL,
        typename MatrixVR, typename VectorS, typename VectorSEP >
inline typename std::enable_if<! detail::is_workspace< VectorSEP >::value,
        std::ptrdiff_t >::type
trsna( const char job, const char howmny, const VectorSELECT& select,
        const MatrixT& t, const MatrixVL& vl, const MatrixVR& vr,
        const VectorS& s, const VectorSEP& sep, const fortran_int_t mm,
        fortran_int_t& m ) {
    return trsna_impl< typename bindings::value_type<
            MatrixT >::type >::invoke( job, howmny, select, t, vl, vr, s, sep,
            mm, m, optimal_workspace() );
}

} // namespace lapack
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
