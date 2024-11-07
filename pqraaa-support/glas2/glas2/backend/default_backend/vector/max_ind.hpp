//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_backend_default_backend_vector_algorithm_max_ind_hpp
#define glas2_backend_default_backend_vector_algorithm_max_ind_hpp

#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/backend/default_backend/default_backend.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cmath>

namespace glas2 {

  template <typename X>
  typename std::enable_if< is<DenseVector,X>::value
                         , typename X::size_type
                         >::type max_ind( default_backend, X const& x ) {
    typename X::value_type max = 0 ;
    typename X::size_type ind = 0 ;
    if (x.size()>0) max = x(0) ;
    for (typename X::size_type i=1; i<x.size(); ++i) {
      if (x(i) > max) {
        max = x(i) ;
        ind = i ;
      }
    }
    return ind ;
  }

} // namespace glas2

#endif
