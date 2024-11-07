//  (C) Copyright Karl Meerbergen & Dries De Samblanx, 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_concept_solve_hpp
#define cork_concept_solve_hpp

#include <type_traits>
#include <cork/vector.hpp>
#include <cork/concept/arithmetic.hpp>

#ifdef CORK_USE_CONCEPTS
#include <concepts>
#endif

namespace CORK {

#ifdef CORK_USE_CONCEPTS
  template <typename F, typename ValueType>
  concept Solve = requires( F solve, ValueType sigma, CORK::vector<ValueType> x, bool is_shift ) {
    requires Arithmetic<ValueType> ;
    {solve(sigma, x, is_shift)} -> std::same_as<void> ;
  } ;

#endif

} // namespace CORK

#endif

