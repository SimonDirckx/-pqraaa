//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_concept_conj_hpp
#define glas2_concept_conj_hpp

#include <type_traits>

#ifndef GLAS_COMPLEX
#define GLAS_COMPLEX
#endif

#ifdef GLAS_COMPLEX
#include <complex>
#endif

namespace glas2 {

  template <typename X>
  typename std::enable_if< std::is_arithmetic<X>::value, X >::type conj( X const& x ) { return x ; }

#ifdef GLAS_COMPLEX
  template <typename X>
  typename std::enable_if< std::is_arithmetic<X>::value,std::complex< X > >::type conj( std::complex<X> const& x ) { return std::conj(x) ; }
#endif

}

#endif
