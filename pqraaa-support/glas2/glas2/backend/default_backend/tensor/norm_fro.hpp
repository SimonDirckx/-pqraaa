//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_backend_default_backend_tensor_algorithm_norm_fro_hpp
#define glas2_backend_default_backend_tensor_algorithm_norm_fro_hpp

#include <glas2/tensor/concept/dense_tensor.hpp>
#include <glas2/backend/default_backend/tensor/iterator.hpp>
#include <glas2/backend/default_backend/default_backend.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cmath>

namespace glas2 {

  template <typename X>
  typename std::enable_if< is<DenseTensor,X>::value && std::is_floating_point<typename X::value_type>::value
                         , decltype( std::abs(typename X::value_type()) )
                         >::type norm_fro( default_backend, X const& x ) {
    decltype( std::abs(typename X::value_type()) ) sum = 0 ;

    for ( tensor_detail::iterator< typename X::shape_type > it(x.shape()); !it.is_end(); ++it ) {
      auto v = x( *it ) ;
      sum += v*v ;
    }
    return std::sqrt( sum ) ;
  }
} // namespace glas2

#endif
