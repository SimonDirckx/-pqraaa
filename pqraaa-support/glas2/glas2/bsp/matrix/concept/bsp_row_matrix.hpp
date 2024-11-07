//  (C) Copyright Karl Meerbergen & Albert-Jan Yzelman 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_bsp_matrix_concept_bsp_row_matrix_hpp
#define glas2_bsp_matrix_concept_bsp_row_matrix_hpp

#include <glas2/bsp/matrix/concept/bsp_matrix.hpp>

namespace glas2 { namespace bsp {

  // Row distributed matrix
  struct BSPRowMatrix
  : BSPMatrix
  {
    typedef BSPRowMatrix type ;
  } ;

} } // namespace glas::bsp

#endif
