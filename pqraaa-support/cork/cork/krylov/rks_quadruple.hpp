//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_krylov_rks_quadruple_hpp
#define cork_krylov_rks_quadruple_hpp

#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <glas2/bindings/vector.hpp>
#include <glas2/bindings/matrix.hpp>
#include <boost/numeric/bindings/lapack/computational/ormqr.hpp>
#include <boost/numeric/bindings/lapack/computational/unmqr.hpp>
#include <boost/numeric/bindings/lapack/computational/geqrf.hpp>
#include <boost/numeric/bindings/herm.hpp>
#include <boost/numeric/bindings/trans.hpp>
#include <type_traits>
#include <complex>

namespace CORK { namespace krylov {

  namespace detail {
    template <typename T, typename EnableIf=void>
    struct rks_quadruple_qr {} ;

    template <typename T>
    struct rks_quadruple_qr<T, typename std::enable_if< std::is_same<T,double>::value || std::is_same<T,float>::value >::type > {
      template <typename A, typename Tau, typename C>
      static int apply_q_right( const A& a, const Tau& tau, C c ) {
        using namespace boost::numeric::bindings ;
        return lapack::ormqr( tag::right(), a, tau, c, lapack::optimal_workspace() ) ;
      }

      template <typename A, typename Tau, typename C>
      static int apply_q_left( const A& a, const Tau& tau, C c ) {
        using namespace boost::numeric::bindings ;
        return lapack::ormqr( tag::left(), trans(a), tau, c, lapack::optimal_workspace() ) ;
      }
    } ;

    template <typename T>
    struct rks_quadruple_qr<std::complex<T>, typename std::enable_if< std::is_same<T,double>::value || std::is_same<T,float>::value >::type > {
      template <typename A, typename Tau, typename C>
      static int apply_q_right( const A& a, const Tau& tau, C c ) {
        using namespace boost::numeric::bindings ;
        return boost::numeric::bindings::lapack::unmqr( tag::right(), a, tau, c, lapack::optimal_workspace() ) ;
      }

      template <typename A, typename Tau, typename C>
      static int apply_q_left( const A& a, const Tau& tau, C c ) {
        using namespace boost::numeric::bindings ;
        return boost::numeric::bindings::lapack::unmqr( tag::left(), herm(a), tau, c, lapack::optimal_workspace() ) ;
      }
    } ;

  }

  template <typename T>
  class rks_quadruple {
    public:
      typedef T                                  value_type ;
      typedef glas2::matrix< value_type >        matrix_type ;
      typedef glas2::vector< value_type >        vector_type ;
      typedef typename matrix_type::size_type    size_type ;

    public:
      rks_quadruple( size_type n, size_type k_max, size_type degree )
      : Q( n, k_max+degree )
      , U( (k_max+degree)*degree, k_max+1 )
      , H( k_max+1, k_max )
      , continuation_combination( k_max, k_max )
      , poles( k_max )
      , rank( 0 )
      , k_max_( k_max )
      , k_( 0 )
      , degree_( degree )
      {
        glas2::fill( U, 0.0 ) ;
        glas2::fill( H, 0.0 ) ;
        glas2::fill( continuation_combination, 0.0 ) ;
        //glas2::fill( Q, 0.0 ) ;
      }

    public:
      glas2::matrix<value_type> K( size_type order ) const {
        glas2::matrix<value_type> K( order+1, order ) ;
        for (size_type i=0; i<order; ++i)
          K( glas2::range(0,order+1), i ) = H( glas2::range(0,order+1), i ) * poles(i) ;
        K(glas2::range(0,order), glas2::all() ) += continuation_combination( glas2::range(0,order), glas2::range(0,order) ) ;
        return std::move(K) ;
      } // K()

      decltype(auto) u_vector( size_type k ) {
        return reshape( U(glas2::all(), k), degree_, k_max_+degree_, glas2::row_major() ) ;
      } // u_vector()

      size_type k() const { return k_ ; }

    public:
      template <typename Alphas, typename Betas>
      void implicit_restart( Alphas const& alphas, Betas const& betas ) {
        assert( alphas.size()==betas.size() ) ;
        assert( alphas.size()<k_ ) ;
        glas2::vector< value_type > tau( k_ ) ;

        glas2::matrxi< value_type > Q_qr( k_+1, k_ ) ;
        Q_qr = glas2::eye(k_1+1,k_) ;

        glas2::matrix<value_type> KH( k_+1, k+ ) ;
        glas2::matrix<value_type> QH( k_+1, 1 ) ;

        for (size_type i=0; i<alphas.size(); ++i) {
          KshiftH( alphas(i), betas(i), KH( glas2::range(0,k_+1), glas2::range(0,k_) ) ) ;

          // Compute Qr factorization
          int info = boost::numeric::bindings::lapack::geqrf( KH, tau ) ;
          assert( info==0 ) ;

          // Apply Q^H to H
          detail::rks_quadruple_qr<value_type>::apply_q_left( KH, tau(glas2::range(0,k_)), H(glas2::range(k_+1,k_)) ) ;

          // Apply Q to Q_qr
          detail::rks_quadruple_qr<value_type>::apply_q_right( KH, tau, Q_qr( glas2::all(), glas2::range(0,k_+1) ) ) ;

          // Make KH upper triangular.

          // Compute Z
          QH( glas2::range(0,k_), 0 ) = H( k_, glas2::range(0,k_) ) ;
          info = boost::numeric::bindings::lapack::geqrf( KH, tau(glas2::range(0,1)) ) ;
          assert( info==0 ) ;

          // Remove first column of H

          // Update poles
          for (size_type j=0; j<k_-1; ++j) poles(j) = poles(j+1) ;
          
          --k_ ;
        }
      } // implicit_restart()

    private:
      // alpha * H + beta * K
      template <typename KH>
      void KshiftH( value_type const& alpha, value_type const& beta, KH K_H ) const {
        assert( K_H.num_rows()==k_+1 ) ;
        assert( K_H.num_columns()==k_ ) ;

        for (size_type i=0; i<k_; ++i)
          K_H( glas2::range(0,k_+1), i ) = H( glas2::range(0,k_+1), i ) * (beta * poles(i) + alpha ) ;
        K_H(glas2::range(0,k_), glas2::all() ) += beta * continuation_combination( glas2::range(0,k_), glas2::range(0,k_) ) ;
      }

    public:
      // step increases by 1.
      void next_step( size_type new_rank ) {
        for (size_type i=rank; i<new_rank; ++i) {
          add_column_of_Q( k_+1 ) ;
        }
        gram_schmidt_u( k_+1 ) ;
        ++k_ ;
      } // next_step()

    private:
      void add_column_of_Q( size_type column ) {
        auto u_v = u_vector( column ) ;

        glas2::shared_vector< value_type > g( rank ) ;
        glas2::range r_rank(0,rank) ;
        g = multiply( transpose(conj(Q(glas2::all(),r_rank))), Q(glas2::all(),rank) ) ;
        Q(glas2::all(),rank) -= multiply( Q(glas2::all(),r_rank), g ) ;
        u_v( glas2::all(), r_rank ) += outer_prod( u_v( glas2::all(), r_rank.size() ), g ) ;

        g = multiply( transpose(conj(Q(glas2::all(),r_rank))), Q(glas2::all(),rank) ) ;
        Q(glas2::all(),rank) -= multiply( Q(glas2::all(),r_rank), g ) ;
        u_v( glas2::all(), r_rank ) += outer_prod( u_v( glas2::all(), r_rank.size() ), g ) ;

        auto norm = norm_2( Q(glas2::all(),rank) ) ;
        if (norm!=0.0) {
          Q(glas2::all(),rank) /= norm ;
          u_v( glas2::all(), rank ) *= norm ;
          ++rank ;
        }

      } // add_column_of_Q()

      void gram_schmidt_u( size_type column ) {
        glas2::vector< value_type > g( column ) ;
        glas2::range r_step(0,column) ;
        g = multiply( transpose(conj(U(glas2::all(),r_step))), U(glas2::all(),column) ) ;
        U(glas2::all(),column) -= multiply( U(glas2::all(),r_step), g ) ;
        H( r_step, column-1 ) = g ;
        g = multiply( transpose(conj(U(glas2::all(),r_step))), U(glas2::all(),column) ) ;
        U(glas2::all(),column) -= multiply( U(glas2::all(),r_step), g ) ;
        H( r_step, column-1 ) += g ;
        H( column, column-1 ) = norm_2( U(glas2::all(),column) ) ;
        if (H(column,column-1)!=0.0) U(glas2::all(),column) /= H(column,column-1) ;
      } // gram_schmidt_u()

    public:
      matrix_type Q ;
      matrix_type U ;
      matrix_type H ;
      matrix_type continuation_combination ;
      vector_type poles ;
      size_type   rank ;

    private:
      size_type k_max_ ;
      size_type k_ ;
      size_type degree_ ;
  } ; // class rks_quadruple
   
} } // namespace CORK::krylov


#endif
