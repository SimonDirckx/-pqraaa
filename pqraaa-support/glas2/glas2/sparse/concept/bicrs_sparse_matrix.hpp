//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_sparse_concept_bicrs_sparse_matrix_hpp
#define glas2_sparse_concept_bicrs_sparse_matrix_hpp

#include <glas2/sparse/concept/sparse_matrix.hpp>

namespace glas2 {

  struct BiCRSSparseMatrix
  : SparseMatrix
  {
    typedef BiCRSSparseMatrix type ;
  } ;

  // concept SparseMatrix
  // structure (Todo)

  // value_type operator()(size_type i, size_type j)
  // size_type nnz()
  // size_type num_rows()
  // size_type num_columns()
  // static const index_base

  // compressed_rows() -> DenseVector
  // columns() -> DenseVector
  // data() -> DenseVector
  // insert(i,j,v)

  // row(i) -> size_type ??
  // column(i) -> size_type ??

} // namespace glas

#endif
