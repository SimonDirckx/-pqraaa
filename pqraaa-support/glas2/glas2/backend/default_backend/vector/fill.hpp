//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_backend_default_backend_vector_algorithm_fill_hpp
#define glas2_backend_default_backend_vector_algorithm_fill_hpp

#include <glas2/backend/default_backend/default_backend.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>

namespace glas2 {

  template <typename V, typename T>
  typename std::enable_if< is<DenseVector,V>::value, V >::type fill( default_backend, V&& v, T const& value ) {
    for (typename std::decay<V>::type::size_type i=0; i<v.size(); ++i) { v(i) = value ; }
    return v ;
  }

} // namespace glas2

#endif
