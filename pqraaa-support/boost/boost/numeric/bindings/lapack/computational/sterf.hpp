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

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_COMPUTATIONAL_STERF_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_COMPUTATIONAL_STERF_HPP

#include <cassert>
#include <boost/numeric/bindings/begin.hpp>
#include <boost/numeric/bindings/is_mutable.hpp>
#include <boost/numeric/bindings/remove_imaginary.hpp>
#include <boost/numeric/bindings/size.hpp>
#include <boost/numeric/bindings/stride.hpp>
#include <boost/numeric/bindings/value_type.hpp>
#include <type_traits>

//
// The LAPACK-backend for sterf is the netlib-compatible backend.
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
inline std::ptrdiff_t sterf( const fortran_int_t n, float* d, float* e ) {
    fortran_int_t info(0);
    LAPACK_SSTERF( &n, d, e, &info );
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * double value-type.
//
inline std::ptrdiff_t sterf( const fortran_int_t n, double* d, double* e ) {
    fortran_int_t info(0);
    LAPACK_DSTERF( &n, d, e, &info );
    return info;
}

} // namespace detail

//
// Value-type based template class. Use this class if you need a type
// for dispatching to sterf.
//
template< typename Value >
struct sterf_impl {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename VectorD, typename VectorE >
    static std::ptrdiff_t invoke( const fortran_int_t n, VectorD& d,
            VectorE& e ) {
        namespace bindings = ::boost::numeric::bindings;
        static_assert( (std::is_same< typename std::remove_const<
                typename bindings::value_type< VectorD >::type >::type,
                typename std::remove_const< typename bindings::value_type<
                VectorE >::type >::type >::value) );
        static_assert( (bindings::is_mutable< VectorD >::value) );
        static_assert( (bindings::is_mutable< VectorE >::value) );
        assert( bindings::size(e) >= n-1 );
        assert( n >= 0 );
        return detail::sterf( n, bindings::begin_value(d),
                bindings::begin_value(e) );
    }

};


//
// Functions for direct use. These functions are overloaded for temporaries,
// so that wrapped types can still be passed and used for write-access. In
// addition, if applicable, they are overloaded for user-defined workspaces.
// Calls to these functions are passed to the sterf_impl classes. In the 
// documentation, most overloads are collapsed to avoid a large number of
// prototypes which are very similar.
//

//
// Overloaded function for sterf. Its overload differs for
// * VectorD&
// * VectorE&
//
template< typename VectorD, typename VectorE >
inline std::ptrdiff_t sterf( const fortran_int_t n, VectorD& d,
        VectorE& e ) {
    return sterf_impl< typename bindings::value_type<
            VectorD >::type >::invoke( n, d, e );
}

//
// Overloaded function for sterf. Its overload differs for
// * const VectorD&
// * VectorE&
//
template< typename VectorD, typename VectorE >
inline std::ptrdiff_t sterf( const fortran_int_t n, const VectorD& d,
        VectorE& e ) {
    return sterf_impl< typename bindings::value_type<
            VectorD >::type >::invoke( n, d, e );
}

//
// Overloaded function for sterf. Its overload differs for
// * VectorD&
// * const VectorE&
//
template< typename VectorD, typename VectorE >
inline std::ptrdiff_t sterf( const fortran_int_t n, VectorD& d,
        const VectorE& e ) {
    return sterf_impl< typename bindings::value_type<
            VectorD >::type >::invoke( n, d, e );
}

//
// Overloaded function for sterf. Its overload differs for
// * const VectorD&
// * const VectorE&
//
template< typename VectorD, typename VectorE >
inline std::ptrdiff_t sterf( const fortran_int_t n, const VectorD& d,
        const VectorE& e ) {
    return sterf_impl< typename bindings::value_type<
            VectorD >::type >::invoke( n, d, e );
}

} // namespace lapack
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
