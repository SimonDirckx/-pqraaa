//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_concept_is_hpp
#define glas2_concept_is_hpp

#include <glas2/concept/concept.hpp>
#include <type_traits>

namespace glas2 {

  template <typename Concept, typename Class, typename EnableIf=void>
  struct is
  //: boost::is_base_of< traits<Class>::concept, Concept >
  : std::is_base_of< Concept, typename glas_concept< typename std::decay<Class>::type >::type >
  {} ;

} // namespace glas

#endif
