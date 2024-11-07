//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_MAX_SOLVER_GROWTH_FACTOR_HPP
#define CORK_OPTIONS_MAX_SOLVER_GROWTH_FACTOR_HPP

#include<cork/options/krylov_breakdown_tolerance.hpp>
#include<cork/options/value_of.hpp>
#include<limits>
#include<cmath>
#include<tuple>

namespace CORK { namespace options {

  template <typename T>
  class max_solver_growth_factor {
    public:
      max_solver_growth_factor()
      : value_( T(1)/std::sqrt(std::numeric_limits<T>::epsilon()) )
      {}

      template <typename ...Ts>
      max_solver_growth_factor( std::tuple<Ts...> const& options )
      : value_( T(1)/ value_of<krylov_breakdown_tolerance<T>>(options) )
      {}

      max_solver_growth_factor( int v )
      : value_(v)
      {}

      T value() const { return value_ ; }
      T& value() { return value_ ; }

    private:
      T value_ ;
  } ;

} } // CORK::options

#endif
