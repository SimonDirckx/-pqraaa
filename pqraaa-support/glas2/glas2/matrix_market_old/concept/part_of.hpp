//  (C) Copyright Karl Meerbergen (2009).
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_market_concept_part_of_hpp
#define glas2_matrix_market_concept_part_of_hpp

#ifdef GLAS_COMPLEX
#include <glas2/valuetype/complex.hpp>
#endif
#include <cassert>

namespace glas2 { namespace matrix_market { 

template <typename X, typename EnableIf=void>
struct part_of_functor {
  typedef bool result_type ;

  inline result_type operator() ( X const& x, int row, int col ) const { return true ; }
} ; // struct part_of_functor

template <typename X>
bool part_of( X const& x, int row, int column ) {
  return part_of_functor<X>() ( x, row, column ) ;
} // part_of()

} } // namespace glas2::matrix_market

#endif
