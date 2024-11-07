//  (C) Copyright Karl Meerbergen & Dries De Samblanx, 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_vector
#define cork_vector

#include <glas2/vector.hpp>
#include <type_traits>

#include <cork/concept/config.hpp>

#ifdef CORK_USE_CONCEPTS
#include <concepts>
#endif

namespace CORK {

  template <typename T>
  struct is_vector
  : std::false_type
  {} ;

#ifdef CORK_USE_CONCEPTS
  template <typename T>
  concept Vector = is_vector<T>::value ;
#endif

  typedef std::ptrdiff_t vector_size_type ;

  template <typename T>
  using vector = glas2::contiguous_vector<T, vector_size_type>;

  template <typename T>
  struct is_vector< vector<T> >
  : std::true_type
  {} ;

} // namespace CORK

#endif

