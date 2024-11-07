//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_concept_is_no_range_or_slice_hpp
#define glas2_vector_concept_is_no_range_or_slice_hpp

#include <glas2/vector/concept/vector.hpp>
#include <glas2/concept/is.hpp>

namespace glas2 {

  template <typename V>
  struct is_no_range_or_slice
  : is<Vector,V>
  {} ;

} // namespace glas

#endif
