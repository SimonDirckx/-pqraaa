//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_concept_contiguous_dense_matrix_hpp
#define glas2_matrix_concept_contiguous_dense_matrix_hpp

#include <glas2/matrix/concept/dense_matrix.hpp>

namespace glas2 {

  struct ContiguousDenseMatrix
  : DenseMatrix
  {
    typedef ContiguousDenseMatrix type ;
  } ;

  // concept: ContiguousDenseMatrix 

} // namespace glas2

#endif
