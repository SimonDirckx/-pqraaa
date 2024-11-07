//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_basis4cork_barycentric_rational_hpp
#define cork_basis4cork_barycentric_rational_hpp

#include <cork/basis/barycentric_rational.hpp>
#include <cork/basis4cork/barycentric_rational_handle.hpp>
#include <cork/basis4cork/basis4cork.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <string>

namespace CORK { namespace basis4CORK {

  template <typename Weights, typename Nodes, typename I>
  class basis4CORK< basis::barycentric_rational<Weights,Nodes,I> >
  {
    public:
      typedef typename std::decay<Weights>::type weights_type ;
      typedef typename std::decay<Nodes>::type   nodes_type ;
      typedef I                                  size_type ;

      template< typename T>
      using value_type_for = typename basis::barycentric_rational<Weights,Nodes,I>::template value_type_for<T> ;

    public:
      explicit basis4CORK( basis::barycentric_rational<Weights,Nodes,I> const& basis )
      : weights_( basis.weights() )
      , nodes_( basis.nodes() )
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
      // Basis4CORK
      template <typename UM, typename ZM, typename Backend=glas2::default_backend>
      void multiply_M( UM const& U, ZM Z, Backend const& backend=glas2::default_backend() ) const {
        assert(U.num_rows()==weights_.size()+1) ;
        auto Z0 = Z(0,glas2::all()) ;
        //Z0 = -U(0,glas2::all()) + U(weights_.size(),glas2::all()) ;
        assign( backend, Z0, -U(0,glas2::all()) + U(weights_.size(),glas2::all()) ) ;
#ifdef CORK_BARYCENTRIC_RATIONAL_WEIGHTED
        for (size_type i=1; i<weights_.size(); ++i) {
          //Z0 += U(i,glas2::all()) ;
          plus_assign( backend, Z0, U(i,glas2::all()) ) ;
          //Z( i, glas2::all() ) = weights_(i-1) * nodes_(i) * U( i+1, glas2::all() ) - weights_(i)*nodes_(i-1) * U( i, glas2::all() ) ;
          auto Zi = Z( i, glas2::all() ) ;
          assign( backend, Zi, weights_(i-1) * nodes_(i) * U( i+1, glas2::all() ) - weights_(i)*nodes_(i-1) * U( i, glas2::all() ) ) ;
        }
#else
        for (size_type i=1; i<weights_.size(); ++i) {
          //Z(0,glas2::all()) += U(i,glas2::all()) ;
          plus_assign( backend, Z0, U(i,glas2::all()) ) ;
          //Z( i, glas2::all() ) = nodes_(i) * U( i+1, glas2::all() ) - nodes_(i-1) * U( i, glas2::all() ) ;
          auto Zi = Z( i, glas2::all() ) ;
          assign( backend, Zi, nodes_(i) * U( i+1, glas2::all() ) - nodes_(i-1) * U( i, glas2::all() ) ) ;
        }
#endif
      } // multiply_M()

      // Basis4CORK
      template <typename UM, typename ZM, typename Backend=glas2::default_backend>
      void multiply_N( UM const& U, ZM Z, Backend const& backend=glas2::default_backend() ) const {
        fill( Z(0,glas2::all()), 0.0 ) ;
#ifdef CORK_BARYCENTRIC_RATIONAL_WEIGHTED
        for (size_type i=1; i<weights_.size(); ++i) {
          //Z( i, glas2::all() ) = weights_(i-1) * U( i+1, glas2::all() ) - weights_(i) * U( i, glas2::all() ) ;
          auto Zi = Z( i, glas2::all() ) ;
          assign( backend, Zi, weights_(i-1) * U( i+1, glas2::all() ) - weights_(i) * U( i, glas2::all() ) ) ;
        }
#else
        //Z( glas2::range(1,Z.num_rows()), glas2::all() ) = U( glas2::range(2,U.num_rows()), glas2::all() ) - U( glas2::range(1,U.num_rows()-1), glas2::all() ) ;
        auto Zrest = Z( glas2::range(1,Z.num_rows()), glas2::all() ) ;
        assign( backend, Zrest, U( glas2::range(2,U.num_rows()), glas2::all() ) - U( glas2::range(1,U.num_rows()-1), glas2::all() ) ) ;
#endif
      } // multiply_N()

    public:
      template <typename AA>
      typename std::enable_if< glas2::is< glas2::DenseMatrix, AA >::value >::type fill_M( AA A ) const {
        assert( A.num_rows()==weights_.size() ) ;
        assert( A.num_columns()==weights_.size()+1 ) ;
        fill(A,0.) ;
#ifdef CORK_BARYCENTRIC_RATIONAL_WEIGHTED
        A( 0, 0 ) = -1.0 ;
        fill( A( 0, glas2::range_from_end(1,0) ), 1.0 ) ;
        diagonal( A, 0 )(glas2::range_from_end(1,0)) = -weights_( glas2::range_from_end(1,0) ) * nodes_( glas2::range_from_end(0,1) ) ;
        diagonal( A, 1 )(glas2::range_from_end(1,0)) = weights_( glas2::range_from_end(0,1) ) * nodes_( glas2::range_from_end(1,0) ) ;
#else
        assert( false ) ; // Not implemented
#endif
      } // fill_M()

      template <typename AA>
      typename std::enable_if< glas2::is< glas2::DenseMatrix, AA >::value >::type fill_N( AA A ) const {
        assert( A.num_rows()==weights_.size() ) ;
        assert( A.num_columns()==weights_.size()+1 ) ;
        fill(A,0.) ;
#ifdef CORK_BARYCENTRIC_RATIONAL_WEIGHTED
        diagonal( A, 0 )(glas2::range_from_end(1,0)) = -weights_( glas2::range_from_end(1,0) ) ;
        diagonal( A, 1 )(glas2::range_from_end(1,0)) = weights_( glas2::range_from_end(0,1) ) ;
#else
        assert( false ) ; // Not implemented
#endif
      } // fill_N()

    public:
      template< typename T>
      using handle_type = barycentric_rational_handle<T, weights_type const&, nodes_type const&, size_type> ;

      template <typename T>
      handle_type<T> handle() const {
        return handle_type<T>( weights_, nodes_ ) ;
      }

    private:
      weights_type const& weights_ ;
      nodes_type const&   nodes_ ;
  } ; // barycentric_rational

} } // namespace CORK::basis4cork

#endif
