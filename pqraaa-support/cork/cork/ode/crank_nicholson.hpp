//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_ode_crank_nicholson_hpp
#define cork_ode_crank_nicholson_hpp

#include <glas2/matrix.hpp>
#include <cmath>
#include <type_traits>

namespace CORK { namespace ode {

  //
  // For solving
  //
  //  B (dx/dt) = A x + f
  //
  // with A and B from the linearization of
  //
  //  P(d/dt) z = b
  //

  // Scheme:
  //   (A - 2/h * B) x_1 = -(A+2/h*B) x_0 - 2*f
  // with f = -e_1 x b
  //
  template <typename ValueType, typename Linearization>
  class crank_nicholson {
    public:
      typedef typename std::decay< Linearization >::type                                linearization_type ;
      typedef ValueType                                                                 value_type ;
      typedef decltype(std::abs(value_type()))                                          real_value_type ;
      typedef typename linearization_type::template solve_handle_type< value_type >     solve_handle_type ;
      typedef typename linearization_type::template multiply_handle_type< value_type >  multiply_handle_type ;

    public:
      crank_nicholson( Linearization& linearization, real_value_type const& time_step, real_value_type const& start_time )
      : linearization_( linearization )
      , solver_handle_( linearization_.template solve_handle<value_type>() )
      , multiply_handle_( linearization_.template multiply_handle<value_type>() )
      , time_( start_time )
      {
        assert( time_step>0.0) ;
        solver_handle_.shift( 2.0 / time_step ) ;
      }

    public:
      template <typename V>
      auto reshape_vector( V v ) const {
        assert( v.size()==linearization_.size()*linearization_.size_of_basis() ) ;
        return reshape( v, linearization_.size_of_basis(), linearization_.size(), glas2::row_major()) ;
      }

    public:
      // y0 and y1 have size of the linearization
      template <typename Y0, typename Y1, typename RHS>
      void time_step( RHS const& rhs, Y0 const& y0, Y1 y1 ) {
        assert( y0.size() == linearization_.size()*linearization_.size_of_basis() ) ;
        assert( y1.size() == linearization_.size()*linearization_.size_of_basis() ) ;
        assert( rhs.size() == linearization_.size() ) ;

        //auto y0_m = reshape(y0, linearization_.size(), linearization_.size_of_basis(), glas2::column_major()) ;
        //auto y1_m = reshape(y1, linearization_.size(), linearization_.size_of_basis(), glas2::column_major()) ;
        auto y0_m = reshape_vector(y0) ;
        auto y1_m = reshape_vector(y1) ;

        multiply_handle_.multiply_B( y0_m, y1_m ) ;
        y1 *= -2.0 * solver_handle_.shift() ;
        y1_m(0,glas2::all()) += value_type(2.0) * rhs ;
        solver_handle_.solve_without_B( y1_m ) ;
        y1 -= y0 ;

        time_ += delta() ;
      } // time_step()

      real_value_type const& time() const { return time_ ; }
      real_value_type delta() const { return glas2::real( 2.0 / solver_handle_.shift() ) ; }

    private:
      Linearization&       linearization_ ;
      solve_handle_type    solver_handle_ ;
      multiply_handle_type multiply_handle_ ;
      real_value_type      time_ ;
  } ; // class crank_nicholson

} } // namespace CORK::ode


#endif
