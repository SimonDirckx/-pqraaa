//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_concept_random_access_container_hpp
#define glas3_concept_random_access_container_hpp

#include <glas3/concept/container.hpp>
#include <glas3/concept/concept.hpp>

#include <array>
#include <vector>
#include <deque>
#include <map>
#include <unordered_map>

namespace glas3 {

  struct RandomAccessContainer
  : Container
  {
    typedef RandomAccessContainer type ;
  } ;

  template <typename... Ts>
  struct concept<std::array<Ts...>>
  : RandomAccessContainer
  {} ;

  template <typename... Ts>
  struct concept<std::deque<Ts...>>
  : RandomAccessContainer
  {} ;

  template <typename... Ts>
  struct concept<std::map<Ts...>>
  : RandomAccessContainer
  {} ;

  template <typename... Ts>
  struct concept<std::unordered_map<Ts...>>
  : RandomAccessContainer
  {} ;

  template <typename... Ts>
  struct concept<std::vector<Ts...>>
  : RandomAccessContainer
  {} ;

} // namespace glas3

#endif
