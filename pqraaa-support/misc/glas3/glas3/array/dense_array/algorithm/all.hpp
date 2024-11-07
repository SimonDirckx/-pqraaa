//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_algorithm_all_hpp
#define glas3_array_dense_array_algorithm_all_hpp

#include <glas3/array/dense_array/concept/dense_array.hpp>
#include <glas3/concept/is.hpp>

#include <type_traits>

namespace glas3 {

  template <typename X>
  typename std::enable_if< is<DenseArray, X>::value, bool >::type
  all( X const& x ) {
    for (typename X::size_type i = 0; i < x.size(); ++i) {
      if ( x[i] == 0 ) return false ;
    }
    return true ;
  }
} // namespace glas3

#endif
