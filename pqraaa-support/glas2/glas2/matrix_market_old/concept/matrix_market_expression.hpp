//  (C) Copyright Karl Meerbergen (2009).
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_market_concept_matrix_market_expression_hpp
#define glas2_matrix_market_concept_matrix_market_expression_hpp

#include <glas2/matrix/concept/matrix.hpp>
#include <type_traits>


namespace glas2 { namespace matrix_market { 

  struct MatrixMarketExpression
  : ::glas2::Matrix
  {
    typedef MatrixMarketExpression type ;
  } ;

} } // namespace glas2::matrix_market

#endif
