//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_algorithm_constant_hpp
#define glas2_matrix_algorithm_constant_hpp

#include <glas2/matrix/type/constant_matrix.hpp>
#include <type_traits>

namespace glas2 {

  template <typename I, typename T>
  typename std::enable_if< std::is_integral<I>::value, constant_matrix<I,T> >::type constant( I n, I m, T const& value ) {
    return constant_matrix<I,T>(n,m,value) ;
  }

} // namespace glas2

#endif
