//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_tensor_concept_strided_dense_tensor_hpp
#define glas2_tensor_concept_strided_dense_tensor_hpp

#include <glas2/tensor/concept/dense_tensor.hpp>

namespace glas2 {

  struct StridedDenseTensor
  : DenseTensor
  {
    typedef StridedDenseTensor type ;
  } ;

  // concept DenseTensor

} // namespace glas

#endif
