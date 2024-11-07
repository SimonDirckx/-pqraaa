//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_algorithm_norm_1_hpp
#define glas2_vector_algorithm_norm_1_hpp

#include <glas2/backend/current_backend.hpp>
#include <glas2/backend/default_backend/vector/norm_1.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cmath>

namespace glas2 {

  template <typename X>
  typename std::enable_if< is<DenseVector,X>::value
                         , decltype( std::abs(typename X::value_type()) )
                         >::type norm_1( X const& x ) {
    return norm_1( current_backend(), x ) ;
  }
} // namespace glas2

#endif
