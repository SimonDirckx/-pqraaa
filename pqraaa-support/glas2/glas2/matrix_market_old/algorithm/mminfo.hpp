//  (C) Copyright Karl Meerbergen 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_market_algorithm_mminfo_hpp
#define glas2_matrix_market_algorithm_mminfo_hpp

#include <glas2/matrix_market/type/mminformation.hpp>
#include <string>
#include <sstream>

namespace glas2 { namespace matrix_market {

  inline mminformation mminfo( std::istream& is ) {
    return mminformation( is ) ;
  }

} } // namespace glas2::matrix_market

#endif
