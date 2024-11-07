//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_krylov_cork_triple_hpp
#define cork_krylov_cork_triple_hpp

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

  template <typename T>
  class cork_triple {
    public:
      typedef T                                  value_type ;
      typedef glas2::matrix< value_type >        matrix_type ;
      typedef glas2::vector< value_type >        vector_type ;
      typedef typename matrix_type::size_type    size_type ;

    private:
      typedef decltype(std::abs(value_type()))   real_type ;

    public:
      cork_triple( size_type n, size_type k_max, size_type degree )
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
      // step increases by 1.
      void next_step( size_type new_rank, bool q_already_orto=true ) {
        if (!q_already_orto) {
          for (size_type i=rank; i<new_rank; ++i) {
            add_column_of_Q( k_+1 ) ;
          }
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
        fill( H( r_step, column-1 ), 0.0 ) ;
        real_type norm_u0( norm_2( U(glas2::all(),column) ) ;
        real_type norm_u( norm_u0 ) ;
        for (int count=0; count<3 && norm_u>0.707*norm_u0; ++count) {
          g = multiply( transpose(conj(U(glas2::all(),r_step))), U(glas2::all(),column) ) ;
          U(glas2::all(),column) -= multiply( U(glas2::all(),r_step), g ) ;
          H( r_step, column-1 ) += g ;
          norm_u0 = norm_u ;
          norm_u = norm_2( U(glas2::all(),column) ) ;
        }
        if (norm_u==0.0) throw std::runtime_error( "CORK:: triple: added vector has norm zero" ) ;
        H( column, column-1 ) = norm_u ;
        if (norm_u!=0.0) U(glas2::all(),column) /= norm_u ;
      } // gram_schmidt_u()

    public:
      template <typename TQ, typename TZ>
      void transform( TQ const& Tq, TZ const& Tz ) {
      } // transform()

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
  } ; // class cork_triple
   
} } // namespace CORK::krylov


#endif
