//  (C) Copyright Karl Meerbergen & Dries De Samblanx, 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_concept_convertible_hpp
#define cork_concept_convertible_hpp

#include <cork/concept/config.hpp>

#include <type_traits>
#include <complex>

#ifdef CORK_USE_CONCEPTS
#include <concepts>
#endif

namespace CORK {

  template <typename T, typename T2>
  struct is_convertible
  : std::is_convertible< T, T2 >
  {} ;

/*  template <typename T, typename S>
  struct is_convertible< T, std::complex<S> >
  : std::true_type
  {} ;

  template <typename T, typename S>
  struct is_convertible< std::complex<T>, S >
  : std::false_type
  {} ;*/

} // namespace CORK

#endif

