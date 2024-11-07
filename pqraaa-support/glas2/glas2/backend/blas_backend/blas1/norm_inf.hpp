//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_backend_blas_backend_blas1_norm_inf_hpp
#define glas2_backend_blas_backend_blas1_norm_inf_hpp

#include <glas2/backend/blas_backend/blas_backend.hpp>
#include <glas2/backend/blas_backend/blas1/is_vector_container.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/bindings/vector/adaptor.hpp>
#include <glas2/concept/is.hpp>
#include <boost/numeric/bindings/blas/level1/iamax.hpp>
#include <type_traits>
#include <cmath>

namespace glas2 {

  template <typename X>
  typename std::enable_if< is_vector_container<X>::value
                         , decltype( std::abs(typename X::value_type()) )
                         >::type norm_inf( blas_backend, X const& x ) {
    int i = boost::numeric::bindings::blas::iamax( x ) ;
    return std::abs( x(i-1) ) ;
  }
} // namespace glas2

#endif
