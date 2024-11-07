//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_concept_matrix_hpp
#define glas2_matrix_concept_matrix_hpp

#include <glas2/concept/expression.hpp>

namespace glas2 {

  struct Matrix
  : Expression
  {
    typedef Matrix type ;
  } ;

  // concept Matrix
  // size_type
  // value_type
  // num_rows()
  // num_columns()

} // namespace glas

#endif
