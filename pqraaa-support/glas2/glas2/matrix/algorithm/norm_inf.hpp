//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_algorithm_norm_inf_hpp
#define glas2_matrix_algorithm_norm_inf_hpp

#include <glas2/backend/default_backend/matrix/norm_inf.hpp>
#include <glas2/backend/current_backend.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cmath>

namespace glas2 {

  template <typename X>
  typename std::enable_if< is<DenseMatrix,X>::value
                         , decltype( std::abs(typename X::value_type()) )
                         >::type norm_inf( X const& x ) {
    return norm_inf( current_backend(), x ) ;
  }
} // namespace glas2

#endif
