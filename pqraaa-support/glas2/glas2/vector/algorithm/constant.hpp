//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_algorithm_constant_hpp
#define glas2_vector_algorithm_constant_hpp

#include <glas2/vector/type/constant_vector.hpp>
#include <type_traits>

namespace glas2 {

  template <typename I, typename T>
  typename std::enable_if< std::is_integral<I>::value, constant_vector<I,T> >::type constant( I n, T const& value ) {
    return constant_vector<I,T>(n,value) ;
  }

} // namespace glas2

#endif
