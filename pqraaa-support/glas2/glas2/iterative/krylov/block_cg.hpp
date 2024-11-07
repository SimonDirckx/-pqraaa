//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_iterative_krylov_block_cg_hpp
#define glas2_iterative_krylov_block_cg_hpp

#include <glas2/iterative/krylov/options.hpp>
#include <glas2/iterative/krylov/tolerance.hpp>
#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 { namespace iterative {

  // r takes the right-hand side. On output, r contains the residual.
  // Initial x is initial solution
  // Without preconditioner
  template <typename X, typename R, typename Op, typename Report>
  typename std::enable_if< glas2::is<X,glas2::DenseMatrix> >::type cg( Op const& op, X x, R r, Report& report, options const& opt ) {
    typedef typename X::value_type             value_type ;
    typedef matrix< value_type >               container_type ;
    typedef matrix< value_type >               dense_type ;
    typedef typename container_type::size_type size_type ;

    container_type  p_c( x.num_rows(), x.num_columns() ) ;
    container_type  q_c( x.num_rows(), x.num_columns() ) ;
    dense_type      H_c( x.num_columns(), x.num_columns() ) ;
    vector< value_type > g_c( x.num_columns() ;
    auto p = p_c.pass_ref() ;
    auto q = q_c.pass_ref() ;
    auto H = H_c.pass_ref() ;
    auto g = g_c.pass_ref() ;

    assert( r.num_columns()==x.num_columns() ) ;
    assert( x.num_rows() == r.num_rows() ) ;

    // Compute residual
    op( x, q ) ;
    r -= q ;
    auto res_norm_0 = norm_fro( r ) ;
    report( res_norm_0, 0 ) ;

    value_type rho2 ;
    value_type rho1 ;
    for ( decltype(opt.max_mat_vec_) it=0 ; it<opt.max_mat_vec_ ; ++it ) {
      rho1 = multiply( transpose(r), r ) ;

      if (it==0) {
        p = r ;
      } else {
        p = r + (rho1/rho2) * p ;
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
} }

#endif
