//  (C) Copyright Karl Meerbergen & Dries De Samblanx, 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_concept_eigenvalue_selector_hpp
#define cork_concept_eigenvalue_selector_hpp

#include <cork/concept/arithmetic.hpp>
#include <cork/concept/config.hpp>

#ifdef CORK_USE_CONCEPTS

#include <cork/vector.hpp>
#include <type_traits>
#include <concepts>

namespace CORK {

  template <typename T, typename Value>
  concept EigenvalueSelector = requires(T const& s, CORK::vector<Value> e, CORK::vector<int> order) {
    requires Arithmetic<Value> ;
    typename T::value_type ;
    requires std::is_convertible< typename T::value_type, Value >::value ;
    {s.n_wanted_max()} -> std::same_as<int> ;
    {s.sort_eigenvalues( e, order )} -> std::same_as<int> ;
  } ;

  template <typename T, typename Value>
  concept EigenvalueSelectorWithShifts = requires(T const& s, CORK::vector<Value> e) {
    requires EigenvalueSelector<T, Value> ;
    {s.shifts( e )} -> std::same_as<void> ;
  } ;

} // namespace CORK

#endif

#endif

