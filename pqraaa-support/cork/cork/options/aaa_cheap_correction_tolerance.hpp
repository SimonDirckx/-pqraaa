//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_AAA_CHEAP_CORRECTION_TOLERANCE_HPP
#define CORK_OPTIONS_AAA_CHEAP_CORRECTION_TOLERANCE_HPP

#include<limits>
#include<tuple>

namespace CORK { namespace options {

  template <typename T>
  class aaa_cheap_correction_tolerance {
    public:
      aaa_cheap_correction_tolerance()
      : value_( 0. )
      {}

      template <typename ...Ts>
      aaa_cheap_correction_tolerance( std::tuple<Ts...> const& options)
      : value_( 0. )
      {}

      aaa_cheap_correction_tolerance( int v )
      : value_(v)
      {}

      T value() const { return value_ ; }
      T& value() { return value_ ; }

    private:
      T value_ ;
  } ;

} } // CORK::options

#endif
