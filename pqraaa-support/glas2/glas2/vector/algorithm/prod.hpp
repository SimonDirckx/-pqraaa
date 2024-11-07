//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_algorithm_prod_hpp
#define glas2_vector_algorithm_prod_hpp

#include <glas2/backend/current_backend.hpp>
#include <glas2/backend/default_backend/vector/prod.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/concept/plus.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cmath>

namespace glas2 {

  template <typename X>
  typename std::enable_if< is<DenseVector,X>::value
                         , typename X::value_type
                         >::type prod( X const& x ) {
    return prod( current_backend(), x ) ;
  }

} // namespace glas2

#endif
