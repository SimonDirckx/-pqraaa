//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_concept_abs_squared_hpp
#define glas2_concept_abs_squared_hpp

#include <glas2/concept/is_arithmetic.hpp>
#include <type_traits>
#include <cmath>

#ifndef GLAS_COMPLEX
#define GLAS_COMPLEX
#endif

#ifdef GLAS_COMPLEX
#include <complex>
#endif

namespace glas2 {

#ifdef GLAS_COMPLEX

  template <typename X>
  typename std::enable_if< glas2::is_arithmetic<X>::value, decltype(std::norm(X())) >::type abs_squared( X const& x ) { return std::norm(x) ; }

#else

  template <typename X>
  typename std::enable_if< glas2::is_arithmetic<X>::value, X >::type abs_squared( X const& x ) { return x*x ; }

#endif

}

#endif
