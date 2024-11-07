//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_algorithm_structure_hpp
#define glas2_vector_algorithm_structure_hpp

#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>

namespace glas2 {

  template <typename V>
  typename std::enable_if< is<DenseVector,V>::value, typename V::size_type >::type structure( V const& v ) {
    return v.size() ;
  }

} // namespace glas2

#endif
