//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_basis4cork_barycentric_rational_strong_hpp
#define cork_basis4cork_barycentric_rational_strong_hpp

#include <cork/basis/barycentric_rational_strong.hpp>
#include <cork/basis4cork/explicit_matrices.hpp>
#include <cork/basis4cork/basis4cork.hpp>
#include <cork/basis4cork/evaluate.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <string>

namespace CORK { namespace basis4CORK {

  template <typename Weights, typename Nodes, typename I>
  class basis4CORK< basis::barycentric_rational_strong<Weights,Nodes,I> >
  : public explicit_matrices< typename basis::barycentric_rational_strong<Weights,Nodes,I>::minimal_value_type >
  {
    private:
      typedef typename basis::barycentric_rational_strong<Weights,Nodes,I>::minimal_value_type value_type ;

    public:
      typedef typename std::decay<Weights>::type weights_type ;
      typedef typename std::decay<Nodes>::type   nodes_type ;
      typedef I                                  size_type ;

      template< typename T>
      using value_type_for = typename basis::barycentric_rational_strong<Weights,Nodes,I>::template value_type_for<T> ;

    public:
      explicit basis4CORK( basis::barycentric_rational_strong<Weights,Nodes,I> const& basis )
      : explicit_matrices< value_type >( basis.num_terms()-1 )
      , weights_( basis.weights() )
      , nodes_( basis.nodes() )
      {
        // Compute M and N
        auto M = this->M() ;
        auto N = this->N() ;
        fill(M, 0.0) ;
        fill(N, 0.0) ;

        diagonal( M, 1 ) = -weights_(glas2::range_from_end(1,0)) * nodes_(glas2::range_from_end(0,1)) ;
        if (weights_.size()>2) diagonal( M, 2 ) = weights_(glas2::range_from_end(0,2)) * nodes_(glas2::range_from_end(1,1)) ;

        diagonal( N, 1 ) = -weights_(glas2::range_from_end(1,0)) ;
        if (weights_.size()>2) diagonal( N, 2 ) = weights_(glas2::range_from_end(0,2)) ;

        M(M.num_rows()-1,glas2::all()) -= glas2::constant_vector<int,value_type>( M.num_columns(), weights_(weights_.size()-2)*nodes_(nodes_.size()-1) ) ;
        N(N.num_rows()-1,glas2::all()) -= glas2::constant_vector<int,value_type>(M.num_columns(), weights_(weights_.size()-2) ) ;
        M(M.num_rows()-1,0) = weights_(weights_.size()-2)*nodes_(nodes_.size()-1) ;
        N(N.num_rows()-1,0) = weights_(weights_.size()-2) ;
      }


    public:
      // Basis4CORK
      int num_terms() const { return nodes_.size()+1 ; }

      typename std::add_lvalue_reference< typename std::add_const<Weights>::type >::type weights() const& {
        return weights_ ;
      }

      typename std::add_lvalue_reference< typename std::add_const<Nodes>::type >::type nodes() const& {
        return nodes_ ;
      }

    public:
      template <typename Coefs>
      void evaluate( Coefs coefs ) const {
        CORK::basis4CORK::evaluate( *this, coefs ) ;
      }

    private:
      weights_type const& weights_ ;
      nodes_type const&   nodes_ ;
  } ; // barycentric_rational

} } // namespace CORK::basis4cork

#endif
