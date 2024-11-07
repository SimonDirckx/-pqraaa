//  (C) Copyright Karl Meerbergen & Dries De Samblanx, 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_concept_set_of_functions_hpp
#define cork_concept_set_of_functions_hpp

#include <type_traits>
#include <cork/vector.hpp>
#include <cork/concept/arithmetic.hpp>

#ifdef CORK_USE_CONCEPTS
#include <concepts>
#endif

namespace CORK {

#ifdef CORK_USE_CONCEPTS
  template <typename F, typename ValueType>
  concept SetOfScalarFunctions = requires ( F fun, ValueType x, CORK::vector<ValueType> v ) {
    requires Arithmetic<ValueType> ;
    {fun(x, v)} ;
  } ;

#endif

} // namespace CORK

#endif

