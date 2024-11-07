//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_concept_contiguous_dense_vector_hpp
#define glas2_vector_concept_contiguous_dense_vector_hpp

#include <glas2/vector/concept/strided_dense_vector.hpp>

namespace glas2 {

  struct ContiguousDenseVector
  : DenseVector
  {
    typedef ContiguousDenseVector type ;
  } ;

  // concept: ContiguousDenseVector 

} // namespace glas

#endif
