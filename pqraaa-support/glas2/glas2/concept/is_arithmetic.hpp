//  (C) Copyright Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_concept_is_arithmetic_hpp
#define glas2_concept_is_arithmetic_hpp

#include <type_traits>

#ifndef GLAS_COMPLEX
#define GLAS_COMPLEX
#endif

#ifdef GLAS_COMPLEX
#include <complex>
#endif

namespace glas2 {

  template <typename T>
  struct is_arithmetic
  : std::is_arithmetic< T >
  {} ;

#ifdef GLAS_COMPLEX
  template <typename T>
  struct is_arithmetic< std::complex<T> >
  : std::is_arithmetic< T >
  {} ;
#endif

} // namespace glas

#endif
