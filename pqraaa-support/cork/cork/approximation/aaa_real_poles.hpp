//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_aaa_real_poles_hpp
#define cork_approximation_aaa_real_poles_hpp

#include <cork/approximation/aaa_approximation.hpp>
#include <cork/basis4cork/barycentric_rational_real_strong.hpp>
#include <cork/basis4cork/barycentric_rational_real.hpp>
#include <cork/lapack/eig.hpp>
#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace CORK { namespace approximation {

  template <typename AAAApproximation>
  decltype(auto) aaa_real_poles( AAAApproximation const& aaa ) {
    typedef typename AAAApproximation::value_type            value_type ;
    typedef std::complex< decltype(std::abs(value_type())) > complex_value_type ;
    typedef decltype(std::abs(value_type()))                 real_value_type ;

    // Set up linearization
    glas2::matrix< real_value_type > M( aaa.n()-1, aaa.n() ) ;
    //glas2::matrix< real_value_type > M( aaa.n(), aaa.n()+1 ) ;
    glas2::matrix< real_value_type > N( M.num_rows(), M.num_columns() ) ;
    glas2::matrix< complex_value_type > X( M.num_rows(), M.num_rows() ) ;
    glas2::matrix< complex_value_type > Y( M.num_rows(), M.num_rows() ) ;

    basis::barycentric_rational_real_strong< decltype(aaa.weights()), decltype(aaa.nodes()) > basis( aaa.weights(), aaa.nodes() ) ;
    //basis::barycentric_rational_real< decltype(aaa.weights()), decltype(aaa.nodes()) > basis( aaa.weights(), aaa.nodes() ) ;
    basis4CORK::basis4CORK< decltype(basis) > basis4cork( basis ) ;

    basis4cork.fill_M( M ) ;
    basis4cork.fill_N( N ) ;
   // std::cout << CORK::matlab(M, "M") << std::endl ;
   // std::cout << CORK::matlab(N, "N") << std::endl ;

    // Solve eigenvalue problem
    glas2::shared_vector< complex_value_type > poles( M.num_rows() ) ;
    auto Mm = M( glas2::all(), glas2::range_from_end(1,0) ) ;
    auto Nm = N( glas2::all(), glas2::range_from_end(1,0) ) ;
    int info = lapack::eig( Mm, Nm, X, poles ) ;
    assert( info>=0 ) ;
    if (info>0) throw std::runtime_error("QZ method failed for AAA_real_poles") ;

    return poles ;
  } // aaa_real_poles()

} } // namespace CORK::approximation

#endif
