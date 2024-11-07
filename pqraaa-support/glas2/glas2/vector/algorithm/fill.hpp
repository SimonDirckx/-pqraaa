//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_algorithm_fill_hpp
#define glas2_vector_algorithm_fill_hpp

#include <glas2/backend/current_backend.hpp>
#include <glas2/backend/default_backend/vector/fill.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>

namespace glas2 {

  template <typename V, typename T>
  typename std::enable_if< is<DenseVector,V>::value, V >::type fill( V&& v, T const& value ) {
    return fill( current_backend(), v, value ) ;
  }

} // namespace glas2

#endif
