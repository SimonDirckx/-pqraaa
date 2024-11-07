//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_concept_dense_vector_hpp
#define glas2_vector_concept_dense_vector_hpp

#include <glas2/vector/concept/vector.hpp>

namespace glas2 {

  struct DenseVector
  : Vector
  {
    typedef DenseVector type ;
  } ;

  // concept DenseVector
  // operator[](size_type i)
  // operator()(size_type i)

} // namespace glas

#endif
