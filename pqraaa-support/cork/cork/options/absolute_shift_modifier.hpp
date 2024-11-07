//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_ABSOLUTE_SHIFT_MODIFIER_HPP
#define CORK_OPTIONS_ABSOLUTE_SHIFT_MODIFIER_HPP

#include<limits>
#include<cmath>
#include<tuple>

namespace CORK { namespace options {

  template <typename T>
  class absolute_shift_modifier {
    public:
      absolute_shift_modifier()
      : value_( 0.0 )
      {}

      template <typename ...Ts>
      absolute_shift_modifier( std::tuple<Ts...> const& options)
      : value_( 0.0 )
      {}

      absolute_shift_modifier( int v )
      : value_(v)
      {}

      T value() const { return value_ ; }
      T& value() { return value_ ; }

    private:
      T value_ ;
  } ;

} } // CORK::options

#endif
