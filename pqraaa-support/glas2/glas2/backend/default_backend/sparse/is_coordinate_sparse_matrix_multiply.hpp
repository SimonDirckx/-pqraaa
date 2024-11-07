//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_backend_default_backend_sparse_is_coordinate_sparse_matrix_multiply_hpp
#define glas2_backend_default_backend_sparse_is_coordinate_sparse_matrix_multiply_hpp

#include <glas2/concept/is.hpp>
#include <type_traits>

namespace glas2 {

  template <typename E>
  struct is_coordinate_sparse_matrix_multiply
  : std::false_type
  {} ;

} // namespace glas2

#endif
