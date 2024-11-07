//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_CHECK_POLES_HPP
#define CORK_OPTIONS_CHECK_POLES_HPP

#include<tuple>

namespace CORK { namespace options {

  class check_poles {
    public:
      check_poles()
      : value_( false )
      {}

      template <typename ...Ts>
      check_poles( std::tuple<Ts...> const& options)
      : value_( false )
      {}

      check_poles( bool v )
      : value_( v )
      {}

      bool value() const { return value_ ; }
      bool& value() { return value_ ; }

    private:
      bool value_ ;
  } ;

} } // CORK::options

#endif
