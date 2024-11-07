//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_algorithm_min_hpp
#define glas2_vector_algorithm_min_hpp

#include <glas2/vector/algorithm/max.hpp>
#include <glas2/vector/algorithm/ops.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>

namespace glas2 {

  template <typename X>
  typename std::enable_if< is<DenseVector,X>::value
                         , typename X::value_type
                         >::type min( X const& x ) {
    return - max( current_backend(), - x ) ;
  }

} // namespace glas2

#endif
