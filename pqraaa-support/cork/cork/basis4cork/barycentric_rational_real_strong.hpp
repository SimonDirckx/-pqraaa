//  (C) Copyright Karl Meerbergen 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_basis4cork_barycentric_rational_real_strong_hpp
#define cork_basis4cork_barycentric_rational_real_strong_hpp

#include <cork/basis/barycentric_rational_real_strong.hpp>
#include <cork/basis4cork/barycentric_rational_real_matrices.hpp>
#include <cork/basis4cork/explicit_matrices.hpp>
#include <cork/basis4cork/basis4cork.hpp>
#include <cork/basis4cork/evaluate.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <string>

namespace CORK { namespace basis4CORK {

  template <typename Weights, typename Nodes, typename I>
  class basis4CORK< basis::barycentric_rational_real_strong<Weights,Nodes,I> >
  : public explicit_matrices< decltype( std::abs(typename basis::barycentric_rational_real_strong<Weights,Nodes,I>::data_value_type()) ) >
  {
    public:
      typedef typename std::decay<Weights>::type weights_type ;
      typedef typename std::decay<Nodes>::type   nodes_type ;
      typedef I                                  size_type ;

    private:
      typedef decltype( std::abs(typename basis::barycentric_rational_real_strong<Weights,Nodes,I>::data_value_type()) ) value_type ;

    public:
      explicit basis4CORK( basis::barycentric_rational_real<Weights,Nodes,I> const& basis )
      : explicit_matrices< value_type >( basis.nodes().size() )
      , weights_( basis.weights() )
      , nodes_( basis.nodes() )
      {
        int size = basis.nodes().size()+1 ;
        glas2::matrix< value_type > M( size-1, size ) ;
        glas2::matrix< value_type > N( size-1, size ) ;
        barycentric_rational_real_matrices( weights_, nodes_, M, N ) ;
//  std::cout << "M = " << nice_print(M) << std::endl ;
//  std::cout << "N = " << nice_print(N) << std::endl ;

        // Eliminate first row and last but first column or last column depending on whether the last support point is real or complex.
        if (glas2::imag(nodes_(nodes_.size()-1))!=0.0) {
          // The two last nodes are complex conjugates.
          this->M() = M( glas2::range_from_end(1,0), glas2::range_from_end(0,1) ) ; this->M()( glas2::all(), size-2 ) = M( glas2::range_from_end(1,0), size-1 ) ;
          this->N() = N( glas2::range_from_end(1,0), glas2::range_from_end(0,1) ) ; this->N()( glas2::all(), size-2 ) = N( glas2::range_from_end(1,0), size-1 ) ;

          //this->M()( glas2::range_from_end(basis.num_terms()-2,0), glas2::all() ) -= outer_prod( M(glas2::range_from_end(basis.num_terms()-1, 0), basis.num_terms()-1 ), M( 0, glas2::range_from_end(0,1) ) ) ;
          //glas2::range_from_end last_three_rows_old( M.num_rows()-3,0 ) ;
          //glas2::range_from_end last_three_rows( M.num_rows()-4,0 ) ;
          //this->M()( last_three_rows, glas2::range(0, size-2) ) -= outer_prod( M(last_three_rows_old, size-2 ), M( 0, glas2::range_from_end(0,2) ) ) ;
          //this->N()( last_three_rows, glas2::range(0, size-2) ) -= outer_prod( N(last_three_rows_old, size-2 ), M( 0, glas2::range_from_end(0,2) ) ) ;
          this->M()( glas2::all(), glas2::range(0, size-2) ) -= outer_prod( M(glas2::range_from_end(1,0), size-2 ), M( 0, glas2::range_from_end(0,2) ) ) ;
          this->N()( glas2::all(), glas2::range(0, size-2) ) -= outer_prod( N(glas2::range_from_end(1,0), size-2 ), M( 0, glas2::range_from_end(0,2) ) ) ;
        } else {
          // Use real barycentric form here.
          this->M() = M( glas2::range_from_end(1,0), glas2::range_from_end(0,1) ) ;
          this->N() = N( glas2::range_from_end(1,0), glas2::range_from_end(0,1) ) ;
          this->M() -= outer_prod( M(glas2::range_from_end(1,0), size-1 ), M( 0, glas2::range_from_end(0,1) ) ) ;
          this->N() -= outer_prod( N(glas2::range_from_end(1,0), size-1 ), M( 0, glas2::range_from_end(0,1) ) ) ;
        }
//  std::cout << "M = " << nice_print(this->M()) << std::endl ;
//  std::cout << "N = " << nice_print(this->N()) << std::endl ;
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
      /*
      template <typename Coefs>
      void evaluate( Coefs coefs ) const {
        CORK::basis4CORK::evaluate( *this, coefs(glas2::range_from_end(nodes_.size(),0) ) ) ;
        coefs(glas2::range_from_end(nodes_.size()+1,0)) *= 
        CORK::basis4CORK::evaluate( *this, coefs(glas2::range(0,nodes_.size()+1)) ) ;
      }*/

    private:
      weights_type const& weights_ ;
      nodes_type const&   nodes_ ;
  } ; // barycentric_rational_real_strong

} } // namespace CORK::basis4cork

#endif
