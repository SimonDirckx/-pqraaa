//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_NUMBER_OF_STAGNATION_RESTARTS_HPP 
#define CORK_OPTIONS_NUMBER_OF_STAGNATION_RESTARTS_HPP

#include<limits>
#include<tuple>

namespace CORK { namespace options {

  class number_of_stagnation_restarts {
    public:
      number_of_stagnation_restarts()
      : value_(-1)
      {}

      template <typename ...Ts>
      number_of_stagnation_restarts( std::tuple<Ts...> const& options)
      : value_(-1)
      {}

      number_of_stagnation_restarts( int v )
      : value_(v)
      {}

      int value() const { return value_ ; }
      int& value() { return value_ ; }

    private:
      int value_ ;
  } ;

} } // CORK::options

#endif
