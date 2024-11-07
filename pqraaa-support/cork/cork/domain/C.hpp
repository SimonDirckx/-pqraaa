//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_domain_C_hpp
#define cork_domain_C_hpp

#include <cmath>
#include <complex>
#include <type_traits>
#include <glas2/vector.hpp>

namespace CORK { namespace domain {

  template <typename T>
  class C {
    public:
      typedef T                               value_type ;
      typedef decltype(std::abs(value_type))  real_type ;

    public:
      auto distance( Point const& point ) const { return 0.0 ; }
  } ;

} } // namespace CORK::domain

#endif
