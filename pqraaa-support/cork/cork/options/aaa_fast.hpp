//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_AAA_FAST_HPP
#define CORK_OPTIONS_AAA_FAST_HPP

#include<limits>
#include<tuple>

namespace CORK { namespace options {

  class aaa_fast {
    public:
      aaa_fast()
      : value_( false )
      {}

      template <typename ...Ts>
      aaa_fast( std::tuple<Ts...> const& options)
      : value_( false )
      {}

      aaa_fast( bool v )
      : value_(v)
      {}

      bool value() const { return value_ ; }
      bool& value() { return value_ ; }

    private:
      bool value_ ;
  } ;

} } // CORK::options

#endif
