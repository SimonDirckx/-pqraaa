//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_iterative_krylov_cg_hpp
#define glas2_iterative_krylov_cg_hpp

#include <glas2/iterative/krylov/options.hpp>
#include <glas2/iterative/krylov/tolerance.hpp>
#include <glas2/vector.hpp>
#include <cassert>

namespace glas2 { namespace iterative {

  // r takes the right-hand side. On output, r contains the residual.
  // Initial x is initial solution
  template <typename X, typename R, typename Prec, typename Op, typename Par, typename Report>
  void cg( Op const& op, Prec const& prec, Par par, X x, R r, Report& report, options const& opt ) {
    typedef typename X::value_type value_type ;
    typedef vector< value_type >   container_type ;

    assert( x.size() == r.size() ) ;

    container_type  p_c( x.size() ) ;
    container_type  q_c( x.size() ) ;
    auto p = p_c.pass_ref() ;
    auto q = q_c.pass_ref() ;
    auto& z = q ; // Save storage

    // Compute residual
    op( x, q ) ;
    r -= q ;
    auto res_norm_0 = norm_2( r ) ;
    report( res_norm_0, 0 ) ;

    value_type rho2 ;
    value_type rho1 ;
    for ( decltype(opt.max_mat_vec_) it=0 ; it<opt.max_mat_vec_ ; ++it ) {
      z = r ; prec( z ) ;
      rho1 = inner_prod( conj(r), z ) ;
      if (it==0) {
        p = z ;
      } else {
        p = z + (rho1/rho2) * p ;
      }
      op( p, q ) ;
      value_type alpha = rho1 / inner_prod( conj(p), q ) ;
      x += alpha * p ;
      r -= alpha * q ;

      auto res_norm = norm_2( r ) ;
      report( res_norm, it ) ;
      if ( res_norm < tolerance( opt, res_norm_0 ) ) return ;

      rho2 = rho1 ;
    }
  } // cg


  // r takes the right-hand side. On output, r contains the residual.
  // Initial x is initial solution
  // Without preconditioner
  template <typename X, typename R, typename Op, typename Report, typename Impl>
  void cg( Op const& op, X x, R r, Report& report, options const& opt, Impl impl=current_backend() ) {
    typedef typename X::value_type value_type ;
    typedef vector< value_type >   container_type ;

    assert( x.size() == r.size() ) ;

    container_type  p_c( x.size() ) ;
    container_type  q_c( x.size() ) ;
    auto p = p_c.pass_ref() ;
    auto q = q_c.pass_ref() ;

    // Compute residual
    op( x, q ) ;
    r -= q ;
    auto res_norm_0 = norm_2( impl, r ) ;
    report( res_norm_0, 0 ) ;

    value_type rho2 ;
    value_type rho1 ;
    for ( decltype(opt.max_mat_vec_) it=0 ; it<opt.max_mat_vec_ ; ++it ) {
      rho1 = inner_prod( impl, conj(r), r ) ;
      if (it==0) {
        p = r ;
      } else {
        p = r + (rho1/rho2) * p ;
      }
      op( p, q ) ;
      value_type alpha = rho1 / inner_prod( impl, conj(p), q ) ;
      x += alpha * p ;
      r -= alpha * q ;

      auto res_norm = norm_2( impl, r ) ;
      report( res_norm, it ) ;
      if ( res_norm < tolerance( opt, res_norm_0 ) ) return ;

      rho2 = rho1 ;
    }
  } // cg
} }

#endif
