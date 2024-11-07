//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_domain_R_hpp
#define cork_domain_R_hpp

#include <cmath>
#include <complex>
#include <type_traits>
#include <glas2/vector.hpp>

namespace CORK { namespace domain {

  template <typename T>
  class R {
    public:
      typedef T                value_type ;
      typedef value_type       real_type ;

      static_assert( std::is_same<value_type,decltype( std::abs(T()) )>::value, "CORK::domain::R<T>: T should be a real arithmetic type" ) ;

    public:
      auto distance( Point const& point ) const {
        return std::abs(glas2::imag(point)) ;
      }
  } ;

} } // namespace CORK::domain

#endif
