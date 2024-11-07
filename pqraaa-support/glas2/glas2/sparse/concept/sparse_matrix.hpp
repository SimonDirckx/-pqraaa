//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_sparse_concept_sparse_matrix_hpp
#define glas2_sparse_concept_sparse_matrix_hpp

#include <glas2/matrix/concept/matrix.hpp>

namespace glas2 {

  struct SparseMatrix
  : Matrix
  {
    typedef SparseMatrix type ;
  } ;

  // concept SparseMatrix
  // structure (Todo)

  // value_type operator()(size_type i, size_type j)
  // index_type num_nz()
  // size_type num_rows()
  // size_type num_columns()
  // size_type row(index_type)
  // size_type column(index_type)

} // namespace glas

#endif
