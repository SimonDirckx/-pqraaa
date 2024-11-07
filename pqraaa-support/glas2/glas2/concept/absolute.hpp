//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_concept_absolute_hpp
#define glas2_concept_absolute_hpp

#include <type_traits>
#include <cmath>

#ifndef GLAS_COMPLEX
#define GLAS_COMPLEX
#endif

#ifdef GLAS_COMPLEX
#include <complex>
#endif

namespace glas2 {

  struct absolute {
    template <typename X>
    decltype( std::abs(X()) ) operator() ( X const& x) const { return std::abs(x) ; }
  } ;

}

#endif
