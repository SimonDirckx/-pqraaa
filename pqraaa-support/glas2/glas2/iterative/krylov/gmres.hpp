//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_iterative_krylov_gmres_hpp
#define glas2_iterative_krylov_gmres_hpp

#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <glas2/scalar.hpp>
#include <glas2/iterative/krylov/tolerance.hpp>
#include <glas2/iterative/krylov/options.hpp>
#include <cassert>
#include <algorithm>
#include <cmath>

namespace glas2 { namespace iterative {

  struct gmres_options
  : options
  {
    inline gmres_options()
    : n_vectors_( 20 )
    {}

    unsigned int n_vectors_ ;
  } ;

  // Op and Prec are here separate operators
/*  template <typename X, typename Y, typename Prec, typename Op, typename Report>
  void gmres( Op const& op, Prec const& prec, X& x, Y const& y, Report& report, gmres_options const& opt ) {
    typedef typename X::value_type      value_type ;
    typedef shared_vector< value_type > vector_type ;
    typedef shared_matrix< value_type > basis_type ;

    basis_type v( x.size(), opt.n_vectors_+1 ) ;
    basis_type h( opt.n_vectors_+1, opt.n_vectors_ ) ;

    basis_type givens( 2, opt.n_vectors_ ) ;
    vector_type s( opt.n_vectors_+1 ) ;

    glas2::fill( h, 0.0 ) ;

    assert( x.size()==y.size() ) ;

    diagonal( h, 0 ) = 1.0 ;

    double res_norm_0 ;

    for ( int it=0; it<=opt.max_mat_vec_/opt.n_vectors_; ++it ) {
      typename column_result_type<basis_type>::type v_0( column( v, 0 ) ) ;
      op( x, v_0 ) ;
      v_0 -= y ;
      prec( v_0 ) ;

      s(0) = norm_2( v_0 ) ;
      v_0 /= s(0) ;
      s(0) = -s(0) ;
      std::fill( begin(s)+1, end(s), 0.0 ) ;

      if (it==0) {
        res_norm_0 = std::abs(s(0)) ;
        report( res_norm_0, 0 ) ;
        if ( res_norm_0== 0.0 ) return ;
      }

      for (unsigned int i=0; i<opt.n_vectors_; ++i ) {
        auto v_i1( v( glas::all(),  i+1 ) ) ;
        op( v(glas::all(),i), v_i1 ) ;
        prec( v_i1 ) ;

        // Perform two loops of classical Gram Schmidt
        {
          vector_type g( i+1 ) ;
          g = multiply( trans( conj( v(glas2::all(), glas2::range(0,i+1) ) ) ), v_i1 ) ;
          h( range(0,i+1), i ) = g ;
          v_i1 -= column_range( v, range(0,i+1) ) * g ;

          g = multiply( transpose( conj( v( glas::all(), range(0,i+1) ) ) ), v_i1 ) ;
          h( range(0,i+1), i ) += g ;
          v_i1 -= multiply( v( all(), range(0,i+1) ), g ) ;
        }
        h( i+1, i ) = norm_2( v_i1 ) ;
        v_i1 /= h( i+1, i ) ;

        // Apply previous Givens rotations.
        for ( int k=0; k<i; ++k ) {
          value_type temp( h(k,i) ) ;
          h(k,i) = temp * givens(0,k) + h(k+1,i) * givens(1,k) ;
          h(k+1,i) = h(k+1,i) * std::conj( givens(0,k) ) - temp * givens(1,k) ;
        }

        // New Givens rotation
        double tt = norm_2( vector_range( column( h, i ), range(i,i+2) ) ) ;
        givens(0,i) = glas::conj( h( i, i ) ) / tt ;
        givens(1,i) = h( i+1, i ) / tt ; // Is real since h( i+1, i ) is real, so, no conj() needed when s is used.
        h( i, i ) = tt ;
        h( i+1, i ) = 0 ;

        s(i+1)  = (-givens(1,i)) * s(i) ;
        s(i)   *= givens(0,i) ;

        report( glas::abs(s(i+1)), it*opt.n_vectors_ + i+1 ) ;
        if ( glas::abs(s(i+1)) <= tolerance( opt, res_norm_0 ) ) {
          range ii(0,i+1) ;
          vector_range(s,ii) *= inverse( upper( column_range( row_range(h,ii), ii ) ) ) ;
          x += column_range(v,ii) * vector_range(s,ii) ;
          return ;
        }
      }
      range im( 0,opt.n_vectors_) ;
      vector_range( s, im ) *= inverse( upper( row_range(h, im ) ) ) ;
      x += column_range( v, im )  * vector_range( s, im ) ;
    }
  } // gmres
*/

  // Op contains preconditioner and matrix
  template <typename X, typename Y, typename Op, typename Par, typename Report>
  void gmres( Op const& op, Par par, X x, Y const& y, Report& report, gmres_options const& opt ) {
    typedef typename X::value_type  value_type ;
    typedef vector< value_type >    vector_type ;
    typedef matrix< value_type >    basis_type ;

    basis_type v_c( x.size(), opt.n_vectors_+1 ) ; auto v = v_c.pass_ref() ;
    basis_type h_c( opt.n_vectors_+1, opt.n_vectors_ ) ; auto h = h_c.pass_ref() ;

    basis_type givens_c( 2, opt.n_vectors_ ) ; auto givens = givens_c.pass_ref() ;
    vector_type s_c( opt.n_vectors_+1 ) ; auto s = s_c.pass_ref() ;

    glas2::fill(h, 0.0) ;

    assert( x.size()==y.size() ) ;

    glas2::fill( diagonal(h,0), 1.0 ) ;

    auto res_norm_0 = std::abs( value_type() ) ;

    for ( int it=0; it<=opt.max_mat_vec_/opt.n_vectors_; ++it ) {
      auto v_0( v(glas2::all(),0) ) ;
      op( x, v_0 ) ;
      v_0 -= y ;

      s(0) = norm_2( v_0 ) ;
      if (it==0) {
        res_norm_0 = std::abs(s(0)) ;
        report( res_norm_0, 0 ) ;
      }

      if ( std::abs(s(0)) <= tolerance( opt, res_norm_0 ) ) {
        return ;
      }

      v_0 /= s(0) ;
      s(0) = -s(0) ;
      glas2::fill( s( glas2::range(1,s.size()) ), 0.0 ) ;

      for ( int i=0; i<opt.n_vectors_; ++i ) {
        auto v_i1( v(glas2::all(), i+1) ) ;
        op( v(glas2::all(),i), v_i1 ) ;

        // Perform two loops of classical Gram-Schmidt
        {
          vector_type g( i+1 ) ;
          g = multiply( transpose( conj( v( glas2::all(), glas2::range(0,i+1) ) ) ), v_i1 ) ;
          h(glas2::range(0,i+1), i) = g ;
          v_i1 -= multiply( v( glas2::all(), glas2::range(0,i+1) ), g ) ;

          g = multiply( transpose( conj( v( glas2::all(), glas2::range(0,i+1) ) ) ), v_i1 ) ;
          h(glas2::range(0,i+1), i ) += g ;
          v_i1 -= multiply( v( glas2::all(), glas2::range(0,i+1) ), g ) ;
        }
        // This can be replaced by modified Gram-Schmidt
        //{
        //  for ( int j=0; j<=i; ++j ) {
        //  g =
        //    h(j, i) = par.innerprod( conj( v( glas2::all(), j ) ), v_i1 ) ;
        //    v_i1 -= v( glas2::all(), j ) * h(j,i) ;
        //  }
        //}
        h( i+1, i ) = par.norm_2( v_i1 ) ;
        v_i1 /= h( i+1, i ) ;

        // Apply previous Givens rotations.
        for ( int k=0; k<i; ++k ) {
          value_type temp( h(k,i) ) ;
          h(k,i) = temp * givens(0,k) + h(k+1,i) * givens(1,k) ;
          h(k+1,i) = h(k+1,i) * glas2::conj( givens(0,k) ) - temp * glas2::conj(givens(1,k)) ;
        }

        // New Givens rotation
        double tt = norm_2( h( glas2::range(i,i+2), i ) ) ;
        givens(0,i) = glas2::conj( h( i, i ) ) / tt ;
        givens(1,i) = h( i+1, i ) / tt ; // s is real
        h( i, i ) = tt ;
        h( i+1, i ) = 0 ;

        s(i+1)  = (-givens(1,i)) * s(i) ;
        s(i)   *= givens(0,i) ;

        report( std::abs(s(i+1)), it*opt.n_vectors_ + i+1 ) ;
        if ( std::abs(s(i+1)) <= tolerance( opt, res_norm_0 ) ) {
          glas2::range ii(0,i+1) ;
          inplace_solve( upper( h(ii,ii) ), s(ii) ) ;
          x += multiply( v(glas2::all(),ii), s(ii) ) ;
          return ;
        }
      }
      glas2::range im( 0,opt.n_vectors_) ;
      inplace_solve( upper( h(im,im) ), s(im) ) ;
      x += multiply( v(glas2::all(),im), s(im) ) ;
    }
  } // gmres

} }

#endif
