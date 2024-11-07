//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_MAX_KRYLOV_DIMENSION_HPP
#define CORK_OPTIONS_MAX_KRYLOV_DIMENSION_HPP

#include<limits>
#include<tuple>

namespace CORK { namespace options {

  class max_krylov_dimension {
    public:
      max_krylov_dimension()
      : value_(0)
      {}

      template <typename ...Ts>
      max_krylov_dimension( std::tuple<Ts...> const& options )
      : value_(0)
      {}

      max_krylov_dimension( int v )
      : value_(v)
      {}

      int value() const { return value_ ; }
      int& value() { return value_ ; }

    private:
      int value_ ;
  } ;

} } // CORK::options

#endif
