//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_iterative_krylov_bicgstab_hpp
#define glas2_iterative_krylov_bicgstab_hpp

#include <glas2/scalar.hpp>
#include <glas2/vector.hpp>
#include <glas2/iterative/krylov/options.hpp>
#include <glas2/iterative/krylov/tolerance.hpp>
#include <string>
#include <cassert>


namespace glas2 { namespace iterative {

  // Op and Prec are here separate operators
  template <typename X, typename Y, typename Prec, typename Par, typename Op, typename Report>
  void bicgstab( Op const& op, Prec const& prec, Par par, X& x, Y const& y, Report& report, options const& opt ) {
    typedef typename X::value_type  value_type ;
    typedef glas2::vector< value_type >    container_type ;

    assert( x.size() == y.size() ) ;

    container_type r_c( x.size() ) ; auto r = r_c.pass_ref() ;
    container_type r_tilde_c( x.size() ) ; auto r_tilde = r_tilde_c.pass_ref() ;
    container_type p_c( x.size() ) ; auto p = p_c.pass_ref() ;
    container_type p_hat_c( x.size() ) ; auto p_hat = p_hat_c.pass_ref() ;
    container_type s_c( x.size() ) ; auto s = s_c.pass_ref() ;
    container_type s_hat_c( x.size() ) ; auto s_hat = s_hat_c.pass_ref() ;
    container_type t_c( x.size() ) ; auto t = t_c.pass_ref() ;
    container_type v_c( x.size() ) ; auto v = v_c.pass_ref() ;

    op( x, r ) ;
    r = y - r ;
    r_tilde = r ;

    double res_norm_0 = norm_2( r ) ;
    report( res_norm_0, 0 ) ;
    if ( res_norm_0== 0.0 ) return ;

    value_type rho_old ;
    value_type rho, alpha, omega ;
    for ( unsigned int it=0 ; it<opt.max_mat_vec_ ; ++it ) {
      rho_old = rho ;
      rho = par.inner_prod( glas2::conj(r_tilde), r ) ;
      if ( rho==0.0 ) throw std::string("Rho is zero") ;

      if (0==it) {
        p = r ;
      } else {
        value_type beta = (rho/rho_old) * (alpha/omega) ;
        p = r + beta * (p - omega*v ) ;
      }

      p_hat = p ;
      prec( p_hat ) ;
      op( p_hat, v ) ;

      alpha = rho / par.inner_prod( glas2::conj(r_tilde), v ) ;
      s = r - alpha * v ;

      double res_norm = par.norm_2(s) ;
      report( res_norm, it*2+1 ) ;
      if ( res_norm < glas2::iterative::tolerance( opt, res_norm_0 ) ) {
        x += alpha * p_hat ;
        return ;
      }

      s_hat = s ;
      prec( s_hat ) ;
      op( s_hat, t ) ;
      omega = par.inner_prod( glas2::conj(t), s ) / par.inner_prod( glas2::conj(t), t ) ;

      x += alpha * p_hat + omega * s_hat ;

      r = s - omega * t ;
      res_norm = par.norm_2( r ) ;
      report( res_norm, it*2+2 ) ;
      if ( res_norm < glas2::iterative::tolerance( opt, res_norm_0 ) ) return ;

      if (omega==0.0) throw std::string("Omega is zero") ;
    }
  } // bicgstab

  // Op contains preconditioner and linear operator.
  template <typename X, typename Y, typename Op, typename Par, typename Report>
  void bicgstab( Op const& op, Par par, X x, Y const& y, Report& report, options const& opt ) {
    typedef typename X::value_type  value_type ;
    typedef glas2::vector< value_type >    container_type ;

    assert( x.size() == y.size() ) ;

    container_type r_tilde_c( x.size() ) ; auto r_tilde = r_tilde_c.pass_ref() ;
    container_type r_c( x.size() ) ; auto r = r_c.pass_ref() ;
    container_type p_c( x.size() ) ; auto p = p_c.pass_ref() ;
    container_type s_c( x.size() ) ; auto s = s_c.pass_ref() ;
    container_type t_c( x.size() ) ; auto t = t_c.pass_ref() ;
    container_type v_c( x.size() ) ; auto v = v_c.pass_ref() ;

    op( x, r ) ;
    r = y - r ;
    r_tilde = r ;

    double res_norm_0 = par.norm_2( r ) ;
    report( res_norm_0, 0 ) ;
    if ( res_norm_0== 0.0 ) return ;

    value_type rho_old ;
    value_type rho=0.0, alpha=0.0, omega=0.0 ;
    for ( unsigned int it=0 ; it<opt.max_mat_vec_ ; ++it ) {
      rho_old = rho ;
      rho = par.inner_prod( glas2::conj(r_tilde), r ) ;
      if ( rho==0.0 ) throw std::string("Rho is zero") ;

      if (0==it) {
        p = r ;
      } else {
        value_type beta = (rho/rho_old) * (alpha/omega) ;
        p = r + beta * (p - omega*v ) ;
      }

      op( p, v ) ;

      alpha = par.inner_prod( glas2::conj(r_tilde), v ) ;
      if (alpha==0.0) throw std::string("Alpha is zero") ;
      alpha = rho / alpha ;
      s = r - alpha * v ;

      double res_norm = par.norm_2(s) ;
      report( res_norm, it*2+1 ) ;
      if ( res_norm < glas2::iterative::tolerance( opt, res_norm_0 ) ) {
        x += alpha * p ;
        return ;
      }

      op( s, t ) ;
      omega = par.inner_prod( glas2::conj(t), t ) ;
      if (omega==0.0) throw std::string("Omega is zero") ;
      omega = par.inner_prod( glas2::conj(t), s ) / omega ; 

      x += alpha * p + omega * s ;

      r = s - omega * t ;
      res_norm = par.norm_2( r ) ;
      report( res_norm, it*2+2 ) ;
      if ( res_norm < glas2::iterative::tolerance( opt, res_norm_0 ) ) return ;

      if (omega==0.0) throw std::string("Omega is zero") ;
    }
  } // bicgstab()
} }

#endif
