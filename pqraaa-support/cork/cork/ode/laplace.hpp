//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_ode_laplace_hpp
#define cork_ode_laplace_hpp

#include <cork/ode/crank_nicholson.hpp>
#include <cork/vector.hpp>
#include <glas2/matrix.hpp>
#include <cmath>
#include <type_traits>

namespace CORK { namespace ode {

  //
  //
  template <typename RealValueType>
  class laplace {
    public:
      typedef RealValueType            value_type ;
      typedef std::complex<value_type> complex_value_type ;
      typedef int                      integer_type ;

    public:
      laplace( CORK::vector<value_type> const& frequencies, CORK::vector<integer_type> const& dofs)
      : frequencies_( frequencies )
      , dofs_( dofs )
      , accum_( frequencies_.size(), dofs_.size() )
      , two_pi( 8.*std::atan(1.0) )
      {
        fill( accum_, 0.0 ) ;
      }

    public:
      template <typename X, typename Stepper>
      void handle_time_step( X const& x, Stepper const& time_stepper ) {
        for (int i=0; i<frequencies_.size(); ++i) {
          accum_( i, glas2::all() ) += ( std::exp( complex_value_type(0.,-two_pi*frequencies_(i))*time_stepper.time() ) * time_stepper.delta() ) * x( dofs_ ) ;
        }
      } // handle_time_step()

      auto const& transform() const { return accum_ ; }

    private:
      //solve_handle_type                 solver_handle_ ;
      //multiply_handle_type              multiply_handle_ ;
      CORK::vector<value_type>          frequencies_ ;
      CORK::vector<integer_type>        dofs_ ;
      glas2::matrix<complex_value_type> accum_ ;
      value_type const two_pi ;
  } ; // class laplace

} } // namespace CORK::ode


#endif
