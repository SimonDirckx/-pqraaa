//  (C) Copyright Karl Meerbergen 2019.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_potential_theory_real_hpp
#define cork_approximation_potential_theory_real_hpp

#include <glas2/vector.hpp>
#include <cassert>
#include <cmath>
#include <limits>

namespace CORK { namespace approximation {

  // Given: border, basis.poles
  template <typename Sigma, typename Xi, typename Basis>
  void potential_theory_leja_bagby_real( Sigma const& sigma, Xi const& xi, Basis& basis ) {
    typedef typename Basis::template value_type< typename Sigma::value_type >          value_type ;
    typedef decltype(std::abs(value_type()))    real_type ;

    glas2::vector< real_type > evals_sigma( 2*sigma.size() ) ;
    glas2::vector< real_type > evals_xi( 2*xi.size() ) ;
    auto evals_sigma_up( evals_sigma( glas2::range(0, sigma.size()) ) ) ;
    auto evals_sigma_low( evals_sigma( glas2::range_from_end( sigma.size(),0 ) ) ) ;
    auto evals_xi_up( evals_xi( glas2::range(0, xi.size()) ) ) ;
    auto evals_xi_low( evals_xi( glas2::range_from_end( xi.size(),0 ) ) ) ;

    // First node
    basis.nodes()(0) = sigma(0) ;
    evals_sigma_up = glas2::abs( sigma - basis.nodes()(0) ) ;
    evals_sigma_low = glas2::abs( sigma - std::conj(basis.nodes()(0)) ) ;
    evals_xi_up = glas2::abs( xi - basis.nodes()(0) ) ;
    evals_xi_low = glas2::abs( xi - std::conj(basis.nodes()(0)) ) ;
    real_type nrm = std::max( norm_inf( evals_sigma ), norm_inf( evals_xi ) ) ;
    evals_sigma /= nrm ;
    evals_xi /= nrm ;

    for (int i=1; i<basis.nodes().size(); ++i) {
      //std::cout << " i = " << i << " : " << evals_sigma(glas2::range(0,std::min<int>(10,evals_sigma.size()))) << " , " << evals_xi(glas2::range(0,std::min<int>(10,evals_sigma.size()))) << std::endl ;
      int sigma_i = glas2::max_ind(evals_sigma) % sigma.size() ;
      int xi_i = glas2::max_ind(evals_xi) % xi.size() ;

      basis.nodes()(i) = sigma( sigma_i ) ;
      basis.poles()(i-1) = xi( xi_i ) ;

//      evals_sigma /= evals_sigma(glas2::max_ind(evals_sigma)) ;
//      evals_xi /= evals_xi(glas2::max_ind(evals_xi)) ;

      if (std::imag(basis.nodes()(i))==0.) {
        evals_sigma_up = evals_sigma_up * glas2::abs( sigma - basis.nodes()(i) ) ;
        evals_sigma_low = evals_sigma_low * glas2::abs( sigma - glas2::conj(basis.nodes()(i)) ) ;
        assert( prod( sigma - basis.poles()(i-1) )!=0.0);
        assert( prod( sigma - std::conj(basis.poles()(i-1)) )!=0.0);
        evals_sigma_up /= glas2::abs( sigma - basis.poles()(i-1) ) ;
        evals_sigma_low /= glas2::abs( sigma - std::conj(basis.poles()(i-1)) ) ;

        evals_xi_up = evals_xi_up * glas2::abs( xi - basis.nodes()(i) ) ;
        evals_xi_low = evals_xi_low * glas2::abs( xi - std::conj(basis.nodes()(i)) ) ;
        for (int j=0; j<xi.size(); ++j) {
          real_type t = std::abs( xi(j) - basis.poles()(i-1) ) ;
          if (t==0.) evals_xi(j) = std::numeric_limits<real_type>::infinity() ;
          else evals_xi(j) /= t ;
        }
        for (int j=0; j<xi.size(); ++j) {
          real_type t = std::abs( xi(j) - std::conj(basis.poles()(i-1)) ) ;
          if (t==0.) evals_xi(xi.size()+j) = std::numeric_limits<real_type>::infinity() ;
          else evals_xi(xi.size()+j) /= t ;
        }
      } else {
        evals_sigma_up = evals_sigma_up * glas2::abs_squared( sigma - basis.nodes()(i) ) ;
        evals_sigma_low = evals_sigma_low * glas2::abs_squared( sigma - glas2::conj(basis.nodes()(i)) ) ;
        assert( prod( sigma - basis.poles()(i-1) )!=0.0);
        assert( prod( sigma - std::conj(basis.poles()(i-1)) )!=0.0);
        evals_sigma_up /= glas2::abs_squared( sigma - basis.poles()(i-1) ) ;
        evals_sigma_low /= glas2::abs_squared( sigma - std::conj(basis.poles()(i-1)) ) ;

        evals_xi_up = evals_xi_up * glas2::abs_squared( xi - basis.nodes()(i) ) ;
        evals_xi_low = evals_xi_low * glas2::abs_squared( xi - std::conj(basis.nodes()(i)) ) ;
        for (int j=0; j<xi.size(); ++j) {
          real_type t = std::norm( xi(j) - basis.poles()(i-1) ) ;
          if (t==0.) evals_xi(j) = std::numeric_limits<real_type>::infinity() ;
          else evals_xi(j) /= t ;
        }
        for (int j=0; j<xi.size(); ++j) {
          real_type t = std::norm( xi(j) - std::conj(basis.poles()(i-1)) ) ;
          if (t==0.) evals_xi(xi.size()+j) = std::numeric_limits<real_type>::infinity() ;
          else evals_xi(xi.size()+j) /= t ;
        }

        basis.nodes()(i+1) = std::conj( basis.nodes()(i) ) ;
        basis.poles()(i) = std::conj( basis.poles()(i-1) ) ;
        ++i ;
      }
    }
    // Change the last pole to infinity.
    //basis.poles()(basis.poles().size()-1) = std::numeric_limits<real_type>::infinity() ;
  } // potential_theory_leja_bagby_real()



} } // namespace CORK::approximation

#endif
