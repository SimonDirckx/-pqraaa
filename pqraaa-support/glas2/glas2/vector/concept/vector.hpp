//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_concept_vector_hpp
#define glas2_vector_concept_vector_hpp

#include <glas2/concept/expression.hpp>

namespace glas2 {

  struct Vector
  : Expression
  {
    typedef Vector type ;
  } ;

  // concept Vector
  // size_type
  // value_type
  // size()

} // namespace glas

#endif
