//  (C) Copyright Karl Meerbergen & Dries De Samblanx, 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_concept_domain_hpp
#define cork_concept_domain_hpp

#include <cork/concept/config.hpp>
#include <cork/concept/vector_of.hpp>

#ifdef CORK_USE_CONCEPTS

#include <cork/vector.hpp>
#include <type_traits>
#include <concepts>
#include <cmath>

namespace CORK {

  template <typename T>
  struct real_part_type {
    typedef decltype( std::abs(T()) ) type ;
  } ;

  template <typename T>
  concept Domain = requires(T const& d, int n, typename T::value_type const& p ) {
    typename T::value_type ;
    {d.discretize_coarse(n)} -> VectorOf<typename T::value_type> ;
    {d.discretize(n)} -> VectorOf<typename T::value_type> ;
    {d.distance(p)} -> std::same_as< typename real_part_type<typename T::value_type>::type > ;
  } ;

} // namespace CORK

#endif

#endif

