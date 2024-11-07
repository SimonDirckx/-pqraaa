//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_DOMAIN_DISTANCD_TOLERANCE_HPP
#define CORK_OPTIONS_DOMAIN_DISTANCD_TOLERANCE_HPP

#include<limits>
#include<tuple>

namespace CORK { namespace options {

  template <typename T>
  class domain_distance_tolerance {
    public:
      domain_distance_tolerance()
      : value_( std::numeric_limits<T>::epsilon()*10. )
      {}

      template <typename ...Ts>
      domain_distance_tolerance( std::tuple<Ts...> const& options)
      : value_( std::numeric_limits<T>::epsilon()*10. )
      {}

      domain_distance_tolerance( T v )
      : value_(v)
      {}

      T value() const { return value_ ; }
      T& value() { return value_ ; }

    private:
      T value_ ;
  } ;

} } // CORK::options

#endif
