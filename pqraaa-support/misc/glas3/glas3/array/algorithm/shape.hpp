//  (C) Copyright Sam Corveleyn 2014.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_algorithm_shape_hpp
#define glas3_array_algorithm_shape_hpp

#include <glas3/array/concept/array.hpp>
#include <glas3/concept/is.hpp>
#include <type_traits>

namespace glas3 {

  template <typename X>
  auto shape( X const& x ) -> typename std::enable_if< is<Array, X>::value, decltype( x.shape().shallow_copy() ) >::type {
    return x.shape().shallow_copy() ;
  }

} // namespace glas3

#endif
