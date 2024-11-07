//  (C) Copyright Karl Meerbergen & Dries De Samblanx, 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_concept_scalar_functions_sequence_hpp
#define cork_concept_scalar_functions_sequence_hpp

#include <type_traits>
#include <cork/vector.hpp>
#include <cork/concept/arithmetic.hpp>
#include <cork/concept/stl_sequence_container.hpp>

#ifdef CORK_USE_CONCEPTS
#include <concepts>
#endif

namespace CORK {

#ifdef CORK_USE_CONCEPTS
  template <typename F, typename A>
  concept ScalarFunction = requires ( F fun, A x ) {
    requires Arithmetic<A> ;
    {fun(x)} -> std::same_as<A> ;
  } ;

  template <typename F, typename A>
  concept ScalarFunctionsSequence = requires ( F sequence ) {
    requires STLSequenceContainer<F> ;
    requires ScalarFunction< typename F::value_type, A > ;
  } ;

#endif

} // namespace CORK

#endif

