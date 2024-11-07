//  (C) Copyright Karl Meerbergen 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_algorithm_flip_hpp
#define glas2_vector_algorithm_flip_hpp

#include <glas2/concept/is.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/vector/type/flip_type.hpp>

namespace glas2 {

  template <typename X>
  typename std::enable_if< is<DenseVector,X>::value, flip_type<X> >::type flip( X x ) {
    return flip_type< X >( x ) ;
  }

} // namespace glas2

#endif
