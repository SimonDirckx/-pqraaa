//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_algorithm_is_equal_hpp
#define glas2_vector_algorithm_is_equal_hpp

#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>

namespace glas2 {

  template <typename X, typename Y>
  typename std::enable_if< is<DenseVector,X>::value && is<DenseVector,Y>::value
                         , bool
                         >::type is_equal( X const& x, Y const& y ) {
    if (x.size()!=y.size()) return false ;
    for (typename X::size_type i=0; i<x.size(); ++i) {
      if (x(i)!=y(i)) return false ;
    }
    return true ;
  }
} // namespace glas2

#endif
