//  (C) Copyright Karl Meerbergen & Dries De Samblanx, 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_concept_vector
#define cork_concept_vector

#include <cork/vector.hpp>
#include <type_traits>

#include <cork/concept/config.hpp>

#ifdef CORK_USE_CONCEPTS
#include <concepts>
#endif

namespace CORK {

#ifdef CORK_USE_CONCEPTS
  template <typename V, typename T>
  concept VectorOf = std::is_base_of<CORK::vector<T>, V>::value ;
#endif

} // namespace CORK

#endif

