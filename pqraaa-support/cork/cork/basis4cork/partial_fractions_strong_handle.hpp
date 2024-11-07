//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_basis4cork_partial_fractions_strong_handle_hpp
#define cork_basis4cork_partial_fractions_strong_handle_hpp

#include <cork/basis4cork/basis4cork.hpp>
#include <cork/exception/linear_solver_failure.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <glas2/bindings/matrix.hpp>
#include <glas2/bindings/vector.hpp>
#include <cassert>
#include <string>
#ifndef NDEBUG
#include <limits>
#endif

namespace CORK { namespace basis4CORK {

  template <typename ValueType, typename Poles, typename Weights, typename I>
  class partial_fractions_strong_handle
  {
    public:
      typedef typename std::decay<Poles>::type   poles_type ;
      typedef typename std::decay<Weights>::type weights_type ;
      typedef I                                  size_type ;

      typedef typename std::common_type< ValueType, typename std::common_type< typename poles_type::value_type, typename weights_type::value_type >::type >::type value_type ;

    private:
 //     typedef glas2::vector< value_type > temp_type ;
      typedef glas2::matrix< value_type > matrix_temp_type ;
      typedef glas2::vector< int > int_vector_type ;
    
    public:
      explicit partial_fractions_strong_handle( Poles poles, Weights weights )
      : poles_( poles )
      , weights_( weights )
      {}

    public:
      // Basis4CORK
      int num_terms() const { return poles_.size()+1 ; }
      int size() const { return num_terms()-1 ; }

    public:
      // Basis4CORK
      void shift( value_type s ) {
        if (prod(poles_-s)==0.) throw exception::linear_solver_failure() ;
        s_ = s ;
      }

      // Basis4CORK
      value_type shift() const { return s_ ; }

      // Basis4CORK
      value_type phi_0() const { return 1. ; }

      // Basis4CORK
      template <typename Z0, typename ZZ>
      void lower_solve_right_hand_side( Z0 const& z0, ZZ Z ) const {
        fill( Z, 0.0 ) ;
        Z( 0, glas2::all() ) = weights_(1)*(s_-poles_(1)) * z0 ;
      } // lower_solve_right_hand_side()

      // Basis4CORK
      template <typename ZM>
      void solve( ZM Z ) const {
        for (typename ZM::size_type i=0; i<Z.num_rows(); ++i) {
          Z( i, glas2::all() ) /= weights_(i)*(poles_(i) - s_) ;
          if (i+1<Z.num_rows()) Z( i+1,glas2::all() ) -= weights_(i+2)*(s_-poles_(i+2)) ;
        }
      } // solve()

    private:
      Poles       poles_ ;
      Weights     weights_ ;
      ValueType   s_ ;
  } ; // partial_fractions_handle

} } // namespace CORK::basis4cork

#endif
