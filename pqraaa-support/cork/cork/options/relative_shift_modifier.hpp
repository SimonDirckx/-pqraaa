//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_RELATIVE_SHIFT_MODIFIER_HPP
#define CORK_OPTIONS_RELATIVE_SHIFT_MODIFIER_HPP

#include<tuple>
#include<limits>
#include<cmath>
#include<cork/options/krylov_breakdown_tolerance.hpp>
#include<cork/options/value_of.hpp>

namespace CORK { namespace options {

  template <typename T>
  class relative_shift_modifier {
    public:
      relative_shift_modifier()
      : value_( std::sqrt( std::sqrt( std::numeric_limits<T>::epsilon() ) ) )
      {}

      template <typename ...Ts>
      relative_shift_modifier( std::tuple<Ts...> const& options)
      : value_( std::pow( value_of<krylov_breakdown_tolerance<decltype(std::abs(T()))>>(options), 0.25 ) )
      {}

      relative_shift_modifier( T const& v )
      : value_(v)
      {}

      T value() const { return value_ ; }
      T& value() { return value_ ; }

    private:
      T value_ ;
  } ;

} } // CORK::options

#endif
