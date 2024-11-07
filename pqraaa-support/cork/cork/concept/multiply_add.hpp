//  (C) Copyright Karl Meerbergen & Dries De Samblanx, 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_concept_multiply_add_hpp
#define cork_concept_multiply_add_hpp

#include <type_traits>
#include <cork/vector.hpp>
#include <cork/concept/arithmetic.hpp>

#ifdef CORK_USE_CONCEPTS
#include <concepts>
#endif

namespace CORK {

#ifdef CORK_USE_CONCEPTS
  template <typename F, typename ValueType>
  concept MultiplyAddAux = requires( F const& multiply_add, int i, ValueType alpha, CORK::vector<ValueType> x, CORK::vector<ValueType> y ) {
    {multiply_add(i, alpha, x, y)} -> std::same_as<void> ;
  } ;

  template <typename F, typename ValueType>
  concept MultiplyAdd = Arithmetic<ValueType> && MultiplyAddAux<F,ValueType> ;

#endif

} // namespace CORK

#endif

