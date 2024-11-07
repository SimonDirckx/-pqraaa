//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_scalar_concept_scalar_hpp
#define glas2_scalar_concept_scalar_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/expression.hpp>
#include <glas2/concept/is_arithmetic.hpp>
#include <type_traits>

namespace glas2 {

  struct Scalar
  : Expression
  {
    typedef Scalar type ;
  } ;

  // concept Scalar

  template <typename T>
  struct is< Scalar, T >
  : glas2::is_arithmetic<T>
  {} ;

} // namespace glas

#endif
