//  (C) Copyright Karl Meerbergen 2023.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_restart_reduction_factor_HPP 
#define CORK_OPTIONS_restart_reduction_factor_HPP

#include<limits>
#include<tuple>

namespace CORK { namespace options {

  class restart_reduction_factor {
    public:
      restart_reduction_factor()
      : value_(70)
      {}

      template <typename ...Ts>
      restart_reduction_factor( std::tuple<Ts...> const& options)
      : value_(70)
      {}

      restart_reduction_factor( int v )
      : value_(v)
      {}

      int value() const { return value_ ; }
      int& value() { return value_ ; }

    private:
      int value_ ;
  } ;

} } // CORK::options

#endif
