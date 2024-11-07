//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_sv_aaa_poles_real_hpp
#define cork_approximation_sv_aaa_poles_real_hpp

#include <cork/basis4cork/barycentric_rational_real.hpp>
#include <cork/approximation/sv_aaa_approximation.hpp>
#include <cork/lapack/eig.hpp>
#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace CORK { namespace approximation {

  template <typename AAAApproximation>
  decltype(auto) sv_aaa_poles_real( AAAApproximation const& aaa ) {
    typedef typename AAAApproximation::value_type   value_type ;
    typedef decltype(std::abs(value_type()))        real_value_type ;
    typedef std::complex< real_value_type >         complex_value_type ;

    // Set up linearization
    glas2::matrix< real_value_type > M1( aaa.n(), aaa.n()+1 ) ;
    glas2::matrix< real_value_type > N1( aaa.n(), aaa.n()+1 ) ;
    glas2::matrix< real_value_type > N_copy( aaa.n(), aaa.n() ) ;
    //glas2::matrix< value_type > M_copy( aaa.n(), aaa.n() ) ;
    glas2::matrix< complex_value_type > X( aaa.n(), aaa.n() ) ;

    /*
    // Sort the poles and weights: real first
    // Now done in SV_AAA_real
    */

    // Get the matrices from the dual linear basis
    auto basis = basis::make_barycentric_rational_real( aaa.weights(), aaa.nodes() ) ;
    CORK::basis4CORK::basis4CORK< decltype(basis) > basis4cork(basis) ;
    basis4cork.fill_M( M1 ) ;
    basis4cork.fill_N( N1 ) ;
    auto M = M1(glas2::all(), glas2::range_from_end(1,0) ) ;
    auto N = N1(glas2::all(), glas2::range_from_end(1,0) ) ;

    // Solve eigenvalue problem
    glas2::vector<complex_value_type> e( aaa.n() ) ;
    N_copy = N ;
    //M_copy = M ;
    //int info = boost::numeric::bindings::lapack::ggev( 'V', 'V', M, N, e, e_beta, Y, X ) ;
    int info = lapack::eig( M, N, X, e ) ;
    assert( info>= 0) ;
    if (info>0) throw std::runtime_error("QZ method failed for SV_AAA 2 PF") ;

    return e ;
  } // sv_aaa_poles_real()

} } // namespace CORK::approximation

#endif
