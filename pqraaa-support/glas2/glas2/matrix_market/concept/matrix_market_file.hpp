//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_market_concept_matrix_market_file_hpp
#define glas2_matrix_market_concept_matrix_market_file_hpp

#include <glas2/sparse/concept/sparse_matrix.hpp>

namespace glas2 {

  struct MatrixMarketFile
  : Matrix
  {
    typedef MatrixMarketFile type ;
  } ;

  // concept SparseMatrix
  // structure (Todo)

  // value_type operator()(size_type i, size_type j)
  // size_type nnz()
  // size_type num_rows()
  // size_type num_columns()

  // begin() -> const_iterator
  // end() -> const_iterator

} // namespace glas

#endif
