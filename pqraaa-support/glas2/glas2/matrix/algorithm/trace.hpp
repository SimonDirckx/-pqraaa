//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_algorithm_trace_hpp
#define glas2_matrix_algorithm_trace_hpp

#include <glas2/backend/default_backend/matrix/trace.hpp>
#include <glas2/backend/current_backend.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cmath>

namespace glas2 {

  template <typename X>
  typename std::enable_if< is<DenseMatrix,X>::value
                         , typename X::value_type
                         >::type trace( X const& x ) {
    return trace( current_backend(), x ) ;
  }
} // namespace glas2

#endif
