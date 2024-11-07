//  (C) Copyright Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_scalar_container_static_const_scalar_hpp
#define glas2_scalar_container_static_const_scalar_hpp

#include <glas2/concept/is_arithmetic.hpp>

namespace glas2 {

  template <typename T, int Value>
  class static_const_scalar
  {
    public:
      typedef T value_type ;

    public:
      static_const_scalar()
      {}

      // Copy
      static_const_scalar( static_const_scalar const& that )
      {}

    public:
      static_const_scalar& operator=( static_const_scalar const& that ) {
        return *this ;
      }

      static_const_scalar& operator=( static_const_scalar&& that ) = delete ;

    public:
      operator T() const { return T(Value) ; }
  } ;


  template <typename T, int N>
  struct is_arithmetic< glas2::static_const_scalar<T,N> >
  : is_arithmetic<T>
  {} ;

} // namespace glas2

#endif
