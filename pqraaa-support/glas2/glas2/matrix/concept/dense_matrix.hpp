//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_concept_dense_matrix_hpp
#define glas2_matrix_concept_dense_matrix_hpp

#include <glas2/matrix/concept/matrix.hpp>

namespace glas2 {

  struct DenseMatrix
  : Matrix
  {
    typedef DenseMatrix type ;
  } ;

  // concept DenseMatrix
  // structure (Todo)
  struct general_matrix {} ;
  struct upper_triangular_matrix {} ;
  struct lower_triangular_matrix {} ;
  struct unit_upper_triangular_matrix {} ;
  struct unit_lower_triangular_matrix {} ;

  // orientation: row_major, columns_major
  struct row_major {} ;
  struct column_major {} ;
  struct any_orientation {} ;

  // operator()(size_type i, size_type j)

} // namespace glas

#endif
