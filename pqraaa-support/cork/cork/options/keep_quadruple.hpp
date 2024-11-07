//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)
#ifndef CORK_OPTIONS_keep_quadruple_HPP
#define CORK_OPTIONS_keep_quadruple_HPP

#include<tuple>
#include<limits>

namespace CORK { namespace options {

  class keep_quadruple {
    public:
      keep_quadruple()
      : value_( true )
      {}

      template <typename ...Ts>
      keep_quadruple( std::tuple<Ts...> const& options)
      : value_( true )
      {}

      keep_quadruple( bool v )
      : value_( v )
      {}

      bool value() const { return value_ ; }
      bool& value() { return value_ ; }

    private:
      bool value_ ;
  } ;

} } // CORK::options

#endif
