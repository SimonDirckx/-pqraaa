//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_basis4cork_barycentric_rational_handle_hpp
#define cork_basis4cork_barycentric_rational_handle_hpp

#include <cork/basis/barycentric_rational.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <glas2/bindings/matrix.hpp>
#include <glas2/bindings/vector.hpp>
#include <boost/numeric/bindings/lapack/computational/getrf.hpp>
#include <boost/numeric/bindings/lapack/computational/getrs.hpp>
#include <cassert>
#include <string>
#include <stdexcept>


namespace CORK { namespace basis4CORK {

  template <typename ValueType, typename Weights, typename Nodes, typename I>
  class barycentric_rational_handle
  {
    public:
      typedef typename std::decay<Weights>::type weights_type ;
      typedef typename std::decay<Nodes>::type   nodes_type ;
      typedef I                                  size_type ;
      typedef ValueType                          value_type ;

    private:
 //     typedef glas2::vector< value_type > temp_type ;
      typedef typename std::common_type< value_type
                                        , typename std::common_type< typename weights_type::value_type
                                                                   , typename nodes_type::value_type
                                                                   >::type
                                        >::type        inner_value_type ;
      typedef glas2::shared_matrix< inner_value_type > matrix_temp_type ;
      typedef glas2::shared_vector< int >              int_vector_type ;
    
    public:
      explicit barycentric_rational_handle( Weights weights, Nodes nodes )
      : weights_( weights )
      , nodes_( nodes )
      , lu_( nodes_.size(), nodes_.size() )
      , pivots_( nodes_.size() )
      {}

    public:
      // Basis4CORK
      int num_terms() const { return nodes_.size()+1 ; }
      int size() const { return num_terms() ; }

      typename std::add_lvalue_reference< typename std::add_const<Weights>::type >::type weights() const& {
        return weights_ ;
      }

      typename std::add_lvalue_reference< typename std::add_const<Nodes>::type >::type nodes() const& {
        return nodes_ ;
      }

    public:
      value_type phi_0() const { return 1. ; }

    public:
      // Basis4CORK
      void shift( value_type s ) {
        s_ = s ;

        fill(lu_,0.0) ;
#ifdef CORK_BARYCENTRIC_RATIONAL_WEIGHTED
        fill( lu_(0,glas2::all()), 1. ) ;
        diagonal( lu_, 0 )(glas2::range_from_end(1,0)) = weights_( glas2::range_from_end(0,1)) * (nodes_(glas2::range_from_end(1,0)) - s_) ;
        diagonal( lu_, -1 ) =  weights_( glas2::range_from_end(1,0)) *(s_ - nodes_(glas2::range_from_end(0,1))) ;
#else
        lu_(0,glas2::all()) = weights_ ;
        diagonal( lu_, 0 )(glas2::range_from_end(1,0)) = nodes_(glas2::range_from_end(1,0)) - s_ ;
        diagonal( lu_, -1 ) = s_ - nodes_(glas2::range_from_end(0,1)) ;
#endif
        int info = boost::numeric::bindings::lapack::getrf( lu_, pivots_ ) ;
        assert( info>=0 ) ;
        if (info!=0) throw std::runtime_error("CORK::barycentric_rational: LU factorization of M-s N failed") ;
      }

      // Basis4CORK
      value_type shift() const { return s_ ; }

      // Basis4CORK
      template <typename Z0, typename ZZ>
      void lower_solve_right_hand_side( Z0 const& z0, ZZ Z ) const {
        assert( Z.num_rows()==num_terms()-1 ) ;
        fill( Z, 0.0 ) ;
        Z( 0, glas2::all() ) = - z0 ;
      } // lower_solve_right_hand_side()

      // Basis4CORK
      template <typename ZM, typename Backend=glas2::default_backend>
      void solve( ZM Z, Backend const& backend=glas2::default_backend() ) const {
        static_assert( std::is_same< inner_value_type, typename ZM::value_type >::value, "CORK::basis4CORK::barycentric_rational_handle: ZM must have the same value_type as the inner_value_type" ) ;

        assert( Z.num_rows()==num_terms()-1 ) ;
        glas2::matrix< typename ZM::value_type > Zc(Z.num_rows(), Z.num_columns() ) ;
        Zc = Z ;
        boost::numeric::bindings::lapack::getrs( lu_, pivots_, Zc ) ;
        Z = Zc ;
      } // solve()

    private:
      Weights             weights_ ;
      Nodes               nodes_ ;
      matrix_temp_type    lu_ ;
      int_vector_type     pivots_ ;
      value_type          s_ ;
  } ; // barycentric_rational

} } // namespace CORK::basis4cork

#endif
