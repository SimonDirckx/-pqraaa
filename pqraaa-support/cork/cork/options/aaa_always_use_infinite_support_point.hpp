//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_aaa_always_use_infinite_support_point_hpp
#define CORK_OPTIONS_aaa_always_use_infinite_support_point_hpp

#include<limits>
#include<tuple>

namespace CORK { namespace options {

  class aaa_always_use_infinite_support_point {
    public:
      aaa_always_use_infinite_support_point()
      : value_( false )
      {}

      template <typename ...Ts>
      aaa_always_use_infinite_support_point( std::tuple<Ts...> const& options)
      : value_( false )
      {}

      aaa_always_use_infinite_support_point( bool v )
      : value_(v)
      {}

      bool const& value() const { return value_ ; }
      bool& value() { return value_ ; }

    private:
      bool value_ ;
  } ;

} } // CORK::options

#endif
