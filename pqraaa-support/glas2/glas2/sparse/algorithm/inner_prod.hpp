//  (C) Copyright Karl Meerbergen 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_sparse_algorithm_inner_prod_hpp
#define glas2_sparse_algorithm_inner_prod_hpp

#include <glas2/backend/current_backend.hpp>
#include <glas2/backend/default_backend/sparse/inner_prod.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/sparse/concept/sparse_vector.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename X, typename Y>
  typename std::enable_if< is<SparseVector,X>::value && is<DenseVector,Y>::value
                         , decltype( typename X::value_type() * typename Y::value_type() )
                         >::type inner_prod( X const& x, Y const& y ) {
    return inner_prod( current_backend(), x, y ) ;
  }

  template <typename X, typename Y>
  typename std::enable_if< is<DenseVector,X>::value && is<SparseVector,Y>::value
                         , decltype( typename X::value_type() * typename Y::value_type() )
                         >::type inner_prod( X const& x, Y const& y ) {
    return inner_prod( current_backend(), y, x ) ;
  }

} // namespace glas2

#endif
