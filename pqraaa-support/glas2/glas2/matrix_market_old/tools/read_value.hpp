//  (C) Copyright Karl Meerbergen 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_toolbox_matrix_market_tools_read_value_hpp
#define glas2_toolbox_matrix_market_tools_read_value_hpp

#include <glas2/matrix_market/type/mminformation.hpp>
#include <string>
#include <sstream>

namespace glas2 { namespace matrix_market {

// valuetype_read
template <typename T>
struct read_value {
  template <typename T2>
  void operator() ( enum field_type const& type, std::istream& s, T2& v ) const {
    if (type==REAL) {
      s >> v ;
    } else if (type==COMPLEX) {
      throw "Complex values cannot be read for this container" ;
    }
  }
} ;

#ifdef GLAS_COMPLEX
template <typename T>
struct read_value< std::complex<T> > {
  template <typename T2>
  void operator() ( enum field_type const& type, std::istream& s, T2& v ) const {
    if (type==REAL) {
      s >> v ;
    } else if (type==COMPLEX) {
      T v1, v2 ;
      s >> v1 >> v2 ;
      v = std::complex<T>( v1, v2 ) ;
    }
  }
} ;
#endif

} } // namespace glas2::matrix_market

#endif
