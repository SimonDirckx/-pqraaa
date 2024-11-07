//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_concept_container_hpp
#define glas3_concept_container_hpp

#include <list>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>

namespace glas3 {

  struct Container
  {
    typedef Container type ;
  } ;

  template <typename... Ts>
  struct concept<std::list<Ts...>>
  : Container
  {} ;

  template <typename... Ts>
  struct concept<std::set<Ts...>>
  : Container
  {} ;

  template <typename... Ts>
  struct concept<std::multiset<Ts...>>
  : Container
  {} ;

  template <typename... Ts>
  struct concept<std::multimap<Ts...>>
  : Container
  {} ;

  template <typename... Ts>
  struct concept<std::unordered_set<Ts...>>
  : Container
  {} ;

  template <typename... Ts>
  struct concept<std::unordered_multiset<Ts...>>
  : Container
  {} ;

  template <typename... Ts>
  struct concept<std::unordered_multimap<Ts...>>
  : Container
  {} ;


} // namespace glas3

#endif
