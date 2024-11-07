//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)
#ifndef CORK_OPTIONS_explicit_projection_HPP
#define CORK_OPTIONS_explicit_projection_HPP

#include<tuple>
#include<limits>

namespace CORK { namespace options {

  class explicit_projection {
    public:
      explicit_projection()
      : value_( false )
      {}

      template <typename ...Ts>
      explicit_projection( std::tuple<Ts...> const& options)
      : value_( false )
      {}

      explicit_projection( bool v )
      : value_( v )
      {}

      bool value() const { return value_ ; }
      bool& value() { return value_ ; }

    private:
      bool value_ ;
  } ;

} } // CORK::options

#endif
