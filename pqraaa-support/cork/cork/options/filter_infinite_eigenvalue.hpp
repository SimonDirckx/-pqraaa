//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_FILTER_INFINITE_EIGENVALUES_HPP
#define CORK_OPTIONS_FILTER_INFINITE_EIGENVALUES_HPP

#include<limits>
#include<tuple>

namespace CORK { namespace options {

  class filter_infinite_eigenvalue {
    public:
      filter_infinite_eigenvalue()
      : value_(0)
      {}

      template <typename ...Ts>
      filter_infinite_eigenvalue( std::tuple<Ts...> const& options )
      : value_(0)
      {}

      filter_infinite_eigenvalue( int v )
      : value_(v)
      {}

      int value() const { return value_ ; }
      int& value() { return value_ ; }

    private:
      int value_ ;
  } ;

} } // CORK::options

#endif
