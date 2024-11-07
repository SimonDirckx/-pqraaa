//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_sv_aaa_poles_hpp
#define cork_approximation_sv_aaa_poles_hpp

#include <cork/approximation/sv_aaa_approximation.hpp>
#include <cork/utility/matlab.hpp>
#include <cork/lapack/eig.hpp>
#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace CORK { namespace approximation {

  struct nofilter {
    template <typename T>
    bool operator() ( T const& t ) const {
      return false ;
    }
  } ;

  template <typename AAAApproximation, typename Filter=nofilter>
  decltype(auto) sv_aaa_poles( AAAApproximation const& aaa, Filter const& filter=nofilter() ) {
    typedef typename AAAApproximation::value_type            value_type ;
    typedef std::complex< decltype(std::abs(value_type())) > complex_value_type ;

    // Set up linearization
    glas2::matrix< complex_value_type > M( aaa.n()+1, aaa.n()+1 ) ;
    glas2::matrix< complex_value_type > N( aaa.n()+1, aaa.n()+1 ) ;
    glas2::matrix< complex_value_type > N_copy( aaa.n()+1, aaa.n()+1 ) ;
    //glas2::matrix< value_type > M_copy( aaa.n(), aaa.n() ) ;

    fill( M, 0.0 ) ;
    fill( N, 0.0 ) ;

    M( 0, glas2::range_from_end(1,0) ) = aaa.weights() ;
    fill( M( glas2::range_from_end(1,0), 0 ), 1.0 ) ;
    diagonal( M, 0 )(glas2::range_from_end(1,0)) = aaa.nodes() ;

    fill( diagonal(N,0)(glas2::range_from_end(1,0)), 1.0 ) ;

    // Solve eigenvalue problem
    glas2::vector<complex_value_type> e( aaa.n()+1 ) ;
    glas2::vector<complex_value_type> e_beta( e.size() ) ;
    N_copy = N ;
    //M_copy = M ;
    int info = boost::numeric::bindings::lapack::ggev( 'N', 'N', M, N, e, e_beta, M, M ) ;
    assert( info>= 0) ;
    if (info>0) throw std::runtime_error("QZ method failed for SV_AAA 2 PF") ;

    // Select wanted finite eigenvalues
    glas2::vector< int > selection_all( e.size()-1 ) ;
    std::cout << e/e_beta << std::endl ;
    int pf_i = glas2::max_ind( glas2::abs(e) / glas2::abs(e_beta) ) ;
    int n_pf = 0 ;
    for (int i=0; i<e.size(); ++i) {
      if (i!=pf_i and !filter(e(i)/e_beta(i))) {
        selection_all( n_pf ) = i ;
        ++n_pf ;
      }
    }
    auto selection( selection_all(glas2::range(0,n_pf)) ) ;

    glas2::shared_vector< complex_value_type > poles( selection.size() ) ;
    poles = e(selection) / e_beta(selection) ;

    return poles ;
  } // sv_aaa_poles()

} } // namespace CORK::approximation

#endif
