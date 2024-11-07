//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_algorithm_assign_hpp
#define glas2_vector_algorithm_assign_hpp

#include <glas2/backend/default_backend/vector/assign.hpp>
#include <glas2/backend/current_backend.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename To, typename From>
  typename std::enable_if< is<DenseVector,To>::value && is<DenseVector,From>::value, To >::type assign( To&& to, From const& from ) {
    return assign( current_backend(), to, from ) ;
  }

} // namespace glas2

#endif
