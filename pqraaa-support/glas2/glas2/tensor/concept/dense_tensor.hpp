//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_tensor_concept_dense_tensor_hpp
#define glas2_tensor_concept_dense_tensor_hpp

#include <glas2/tensor/concept/tensor.hpp>

namespace glas2 {

  struct DenseTensor
  : Tensor
  {
    typedef DenseTensor type ;
  } ;

  // concept Tensor
  // size_type
  // value_type
  // order()
  // sizes()

} // namespace glas

#endif
