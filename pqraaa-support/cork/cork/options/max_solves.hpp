//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_MAX_SOLVES_HPP
#define CORK_OPTIONS_MAX_SOLVES_HPP

#include<limits>
#include<tuple>

namespace CORK { namespace options {

  class max_solves {
    public:
      max_solves()
      : value_(10000)
      {}

      template <typename ...Ts>
      max_solves( std::tuple<Ts...> const& options )
      : value_(10000)
      {}

      max_solves( int v )
      : value_(v)
      {}

      int value() const { return value_ ; }
      int& value() { return value_ ; }

    private:
      int value_ ;
  } ;

} } // CORK::options

#endif
