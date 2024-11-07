//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_AAA_MAX_stagnation_iterations_HPP
#define CORK_OPTIONS_AAA_MAX_stagnation_iterations_HPP

#include<limits>
#include<tuple>

namespace CORK { namespace options {

  class aaa_max_stagnation_iterations {
    public:
      aaa_max_stagnation_iterations()
      : value_(1)
      {}

      template <typename ...Ts>
      aaa_max_stagnation_iterations( std::tuple<Ts...> const& options)
      : value_(1)
      {}

      aaa_max_stagnation_iterations( int v )
      : value_(v)
      {}

      int value() const { return value_ ; }
      int& value() { return value_ ; }

    private:
      int value_ ;
  } ;

} } // CORK::options

#endif
