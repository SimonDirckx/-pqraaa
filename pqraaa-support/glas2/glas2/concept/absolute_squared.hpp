//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_concept_absolute_squared_hpp
#define glas2_concept_absolute_squared_hpp

#include <glas2/concept/abs_squared.hpp>
#include <type_traits>
#include <cmath>

namespace glas2 {

#ifndef GLAS_COMPLEX
#define GLAS_COMPLEX
#endif

#ifdef GLAS_COMPLEX

  struct absolute_squared {
    template <typename X>
    decltype( glas2::abs_squared(X()) ) operator() ( X const& x) const { return glas2::abs_squared(x) ; }
  } ;

#else

  struct absolute_squared {
    template <typename X>
    X operator() ( X const& x) const { return x*x ; }
  } ;

#endif

}

#endif
