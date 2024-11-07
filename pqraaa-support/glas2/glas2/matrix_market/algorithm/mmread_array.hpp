//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_market_algorithm_mmread_hpp
#define glas2_matrix_market_algorithm_mmread_hpp

#include <glas2/matrix_market/type/mm_reader.hpp>

namespace glas2 {

  template <typename T>
  inline matrix_market::mm_reader<T> mmread( std::istream& stream ) {
    return matrix_market::mm_reader<T>( stream ) ;
  }

}

#endif
