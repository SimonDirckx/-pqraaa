//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_concept_imag_hpp
#define glas2_concept_imag_hpp

#include <type_traits>

#ifndef GLAS_COMPLEX
#define GLAS_COMPLEX
#endif

#ifdef GLAS_COMPLEX
#include <complex>
#endif

namespace glas2 {

  template <typename X>
  typename std::enable_if< std::is_arithmetic<X>::value, X >::type imag( X const& x ) { return 0. ; }

#ifdef GLAS_COMPLEX
  template <typename X>
  typename std::enable_if< std::is_arithmetic<X>::value, X >::type imag( std::complex<X> const& x ) { return x.imag() ; }
#endif

}

#endif
