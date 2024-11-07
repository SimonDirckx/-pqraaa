//  (C) Copyright Karl Meerbergen & Dries De Samblanx, 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_basis_domain_hpp
#define cork_basis_domain_hpp

#include <cork/concept/config.hpp>

#ifdef CORK_USE_CONCEPTS

#include <cork/vector.hpp>
#include <type_traits>
#include <concepts>
#include <cmath>

namespace CORK { namespace Concept {

  template <typename T>
  concept Basis = requires(T const& basis ) {
    typename T::size_type ;
    {basis.num_terms()} -> std::same_as<typename T::size_type> ;
  } ;

  template <typename T, typename TT>
  concept BasisFor = Basis<T> && requires(T const& basis, TT arg, CORK::vector<TT> values ) {
    {basis.evaluate(arg, values)} -> std::same_as<void> ;
  } ;

} } // namespace CORK::Concept

#endif

#endif

