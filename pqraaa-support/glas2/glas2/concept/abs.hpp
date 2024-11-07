//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_concept_abs_hpp
#define glas2_concept_abs_hpp

#include <type_traits>

#ifndef GLAS_COMPLEX
#define GLAS_COMPLEX
#endif

#include <cmath>
#ifdef GLAS_COMPLEX
#include <complex>
#endif

namespace glas2 {

  template <typename X>
  typename std::enable_if< std::is_arithmetic<X>::value, decltype( std::abs( X() ) ) >::type abs( X const& x ) { return std::abs(x) ; }

#ifdef GLAS_COMPLEX
  template <typename X>
  typename std::enable_if< std::is_arithmetic<X>::value, X >::type abs( std::complex<X> const& x ) { return std::abs(x) ; }
#endif

}

#endif
