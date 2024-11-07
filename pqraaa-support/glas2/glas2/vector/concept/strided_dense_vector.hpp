//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_concept_strided_dense_vector_hpp
#define glas2_vector_concept_strided_dense_vector_hpp

#include <glas2/vector/concept/dense_vector.hpp>

namespace glas2 {

  struct StridedDenseVector
  : DenseVector
  {
    typedef StridedDenseVector type ;
  } ;

  // concept StridedDenseVector
  // size_type stride()
  // value_type* storage_ptr()
  // value_type const* storage_ptr()

} // namespace glas

#endif
