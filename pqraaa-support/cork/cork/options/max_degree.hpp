//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_MAX_DEGREE_HPP
#define CORK_OPTIONS_MAX_DEGREE_HPP

#include<limits>
#include<tuple>

namespace CORK { namespace options {

  class max_degree {
    public:
      max_degree()
      : value_(100)
      {}

      template <typename ...Ts>
      max_degree( std::tuple<Ts...> const& options)
      : value_(100)
      {}

      max_degree( int v )
      : value_(std::max(3,v))
      {}

      int value() const { return value_ ; }
      int& value() { return value_ ; }

    private:
      int value_ ;
  } ;

} } // CORK::options

#endif
