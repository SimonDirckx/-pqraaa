//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_tensor_concept_tensor_unfolding_matrix_hpp
#define glas2_tensor_concept_tensor_unfolding_matrix_hpp

#include <glas2/matrix/concept/dense_matrix.hpp>

namespace glas2 {

  struct TensorUnfoldingMatrix
  : DenseMatrix
  {
    typedef TensorUnfoldingMatrix type ;
  } ;

  // concept DenseMatrix
  // orientation
  // value_type* storage_ptr()
  // value_type const* storage_ptr()

} // namespace glas

#endif
