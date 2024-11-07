//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_backend_blas_backend_blas1_inner_prod_hpp
#define glas2_backend_blas_backend_blas1_inner_prod_hpp

#include <glas2/backend/blas_backend/blas_backend.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/bindings/vector/adaptor.hpp>
#include <glas2/concept/is.hpp>
#include <boost/numeric/bindings/blas/level1/dot.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename X, typename Y>
  typename std::enable_if< is<DenseVector,X>::value && is<DenseVector,Y>::value
                         , decltype( typename X::value_type() * typename Y::value_type() )
                         >::type inner_prod( blas_backend, X const& x, Y const& y ) {
    assert( x.size()==y.size() ) ;
    return boost::numeric::bindings::blas::dot( x, y ) ;
  }
} // namespace glas2

#endif
