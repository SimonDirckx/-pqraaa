//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_concept_semi_contiguous_dense_matrix_hpp
#define glas2_matrix_concept_semi_contiguous_dense_matrix_hpp

#include <glas2/matrix/concept/dense_matrix.hpp>

namespace glas2 {

  struct SemiContiguousDenseMatrix
  : DenseMatrix
  {
    typedef ContiguousDenseMatrix type ;
  } ;

  // concept: SemiContiguousDenseMatrix 
  // size_type leading_dimension()

} // namespace glas

#endif
