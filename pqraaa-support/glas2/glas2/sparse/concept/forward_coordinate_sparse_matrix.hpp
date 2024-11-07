//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_sparse_concept_forward_coordinate_sparse_matrix_hpp
#define glas2_sparse_concept_forward_coordinate_sparse_matrix_hpp

#include <glas2/sparse/concept/coordinate_sparse_matrix.hpp>

namespace glas2 {

  struct ForwardCoordinateSparseMatrix
  : CoordinateSparseMatrix
  {
    typedef ForwardCoordinateSparseMatrix type ;
  } ;

  // begin()
  // end()
  //
  // iterator:
  //   operator++
  //   operator==
  //   operator!=
  //   row()
  //   column()
  //   value()

} // namespace glas

#endif
