//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_tensor_algorithm_norm_fro_hpp
#define glas2_tensor_algorithm_norm_fro_hpp

#include <glas2/backend/default_backend/tensor/norm_fro.hpp>
#include <glas2/backend/current_backend.hpp>
#include <glas2/tensor/concept/dense_tensor.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cmath>

namespace glas2 {

  template <typename X>
  typename std::enable_if< is<DenseTensor,X>::value
                         , decltype( std::abs(typename X::value_type()) )
                         >::type norm_fro( X const& x ) {
    return norm_fro( current_backend(), x ) ;
  }
} // namespace glas2

#endif
