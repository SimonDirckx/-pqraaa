//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_basis4cork_partial_fractions_real_handle_hpp
#define cork_basis4cork_partial_fractions_real_handle_hpp

#include <cork/basis4cork/solve_2_by_2.hpp>
#include <cork/basis4cork/basis4cork.hpp>
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
  class partial_fractions_real_handle
  {
    public:
      typedef typename std::decay<Poles>::type   poles_type ;
      typedef typename std::decay<Weights>::type weights_type ;
      typedef I                                  size_type ;

      typedef  typename std::common_type< ValueType, decltype( std::abs(typename std::common_type< typename poles_type::value_type, typename weights_type::value_type>::type() ) ) >::type value_type ;

    private:
 //     typedef glas2::vector< value_type > temp_type ;
      typedef glas2::matrix< value_type > matrix_temp_type ;
      typedef glas2::vector< int > int_vector_type ;
    
    public:
      explicit partial_fractions_real_handle( Poles poles, Weights weights )
      : poles_( poles )
      , weights_( weights )
      {}

    public:
      // Basis4CORK
      int num_terms() const { return poles_.size()+1 ; }
      int size() const { return num_terms() ; }

    public:
      // Basis4CORK
      void shift( value_type s ) {
        s_ = s ;
      }

      // Basis4CORK
      value_type shift() const { return s_ ; }

      value_type phi_0() const { return 1.0 ; }

      // Basis4CORK
      template <typename Z0, typename ZZ, typename Backend=glas2::default_backend>
      void lower_solve_right_hand_side( Z0 const& z0, ZZ Z, Backend const& backend=glas2::default_backend() ) const {
        for (typename ZZ::size_type i=0; i<Z.num_rows(); ++i) {
          Z( i, glas2::all() ) = std::real(weights_(i)) * z0 ;
          if (std::imag(weights_(i))!=0.0) {
            ++i ;
            Z( i, glas2::all() ) = std::imag(weights_(i-1)) * z0 ;
          }
        }
      } // lower_solve_right_hand_side()

      // Basis4CORK
      template <typename ZM, typename Backend=glas2::default_backend>
      void solve( ZM Z, Backend const& backend=glas2::default_backend() ) const {
        glas2::matrix<value_type> AB(2,2) ;
        for (typename ZM::size_type i=0; i<Z.num_rows(); ++i)
          if (poles_(i).imag()==0.0) {
            Z( i, glas2::all() ) /= (poles_(i).real() - s_) ;
          } else {
            AB(0,0) = poles_(i).real() - s_ ; AB(0,1) = -poles_(i).imag() ;
            AB(1,0) = poles_(i).imag() ; AB(1,1) = poles_(i).real() - s_ ;
            solve_2_by_2( AB, Z(glas2::range(i,i+2), glas2::all()) ) ;
            ++i ;
          }
      } // solve()

    private:
      Poles       poles_ ;
      Weights     weights_ ;
      value_type  s_ ;
  } ; // partial_fractions_real_handle

} } // namespace CORK::basis4cork

#endif
