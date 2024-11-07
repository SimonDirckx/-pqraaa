//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_aaa_poles_hpp
#define cork_approximation_aaa_poles_hpp

#include <cork/approximation/aaa_approximation.hpp>
#include <cork/lapack/eig.hpp>
#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace CORK { namespace approximation {

  template <typename AAAApproximation>
  decltype(auto) aaa_poles( AAAApproximation const& aaa ) {
    typedef typename AAAApproximation::value_type            value_type ;
    typedef std::complex< decltype(std::abs(value_type())) > complex_value_type ;

    glas2::range range_n( 0, aaa.n() ) ;

    // Set up linearization
    glas2::matrix< complex_value_type > M( aaa.n(), aaa.n() ) ;
    glas2::matrix< complex_value_type > N( aaa.n(), aaa.n() ) ;
    glas2::matrix< complex_value_type > N_copy( aaa.n(), aaa.n() ) ;
    //glas2::matrix< value_type > M_copy( aaa.n(), aaa.n() ) ;
    glas2::matrix< complex_value_type > X( aaa.n(), aaa.n() ) ;
    glas2::matrix< complex_value_type > Y( aaa.n(), aaa.n() ) ;

    fill( M, 0.0 ) ;
    fill( N, 0.0 ) ;

    M( 0, glas2::all() ) = aaa.weights() ;
    diagonal( M, -1 ) = - aaa.nodes()(glas2::range_from_end(0,1)) ;
    diagonal( M, 0 )(glas2::range_from_end(1,0)) = aaa.nodes()(glas2::range_from_end(1,0)) ;

    fill( diagonal(N,-1), -1.0 ) ;
    fill( diagonal(N,0)(glas2::range_from_end(1,0)), 1.0 ) ;

    // Solve eigenvalue problem
    glas2::vector<complex_value_type> e( aaa.n() ) ;
    glas2::vector<complex_value_type> e_beta( e.size() ) ;
    N_copy = N ;
    //M_copy = M ;
    int info = boost::numeric::bindings::lapack::ggev( 'V', 'V', M, N, e, e_beta, Y, X ) ;
    assert( info>= 0) ;
    if (info>0) throw std::runtime_error("QZ method failed for SV_AAA 2 PF") ;

/*    for (int i=0; i<e.size(); ++i) {
      std::cout << "Right residual " << norm_2( e_beta(i) * multiply(M_copy, X(glas2::all(),i)) - e(i) * multiply(N_copy, X(glas2::all(),i)) ) << " " 
                << "Left residual " << norm_2( e_beta(i) * multiply(transpose(M_copy), conj(Y(glas2::all(),i))) - e(i) * multiply(transpose(N_copy), conj(Y(glas2::all(),i))) ) << std::endl ;
    }*/

//    std::cout << "Poles " << e / e_beta << std::endl ;

    glas2::shared_vector< complex_value_type > poles( e.size() ) ;
    poles = e/ e_beta ;

    return poles ;
  } // aaa_poles()

} } // namespace CORK::approximation

#endif
