//  (C) Copyright Karl Meerbergen 2009.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_market_algorithm_mmwrite_hpp
#define glas2_matrix_market_algorithm_mmwrite_hpp

#include <glas2/matrix_market/expression/mmwrite_expression.hpp>

#undef GLAS_INCLUDE_LEVEL
#define GLAS_INCLUDE_LEVEL GLAS_INCLUDE_TOOLBOX

namespace glas2 { namespace matrix_market {

  template <typename E>
  inline mmwrite_expression<E> mmwrite( E const& e, std::string const& s = "" ) {
    return mmwrite_expression<E>( e, s ) ;
  }

} } // namespace glas2::matrix_market


#endif
