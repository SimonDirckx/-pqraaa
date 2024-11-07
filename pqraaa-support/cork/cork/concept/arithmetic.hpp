//  (C) Copyright Karl Meerbergen & Dries De Samblanx, 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_concept_arithmetic_hpp
#define cork_concept_arithmetic_hpp

#include <cork/concept/config.hpp>

#include <type_traits>
#include <complex>

#ifdef CORK_USE_CONCEPTS
#include <concepts>
#endif

namespace CORK {

  template <typename T>
  struct is_arithmetic
  : std::false_type
  {} ;

  template<>
  struct is_arithmetic<float>
  : std::true_type
  {} ;

  template<>
  struct is_arithmetic<double>
  : std::true_type
  {} ;

  template <typename T>
  struct is_arithmetic< std::complex<T> >
  : std::is_arithmetic<T>
  {} ;

#ifdef CORK_USE_CONCEPTS
  template <typename T>
  concept Arithmetic = CORK::is_arithmetic<T>::value ;
#endif

} // namespace CORK

#endif

