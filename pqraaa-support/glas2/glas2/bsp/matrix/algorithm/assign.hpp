//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_bsp_matrix_algorithm_assign_hpp
#define glas2_bsp_matrix_algorithm_assign_hpp

#include <glas2/bsp/backend/default/matrix/assign.hpp>
#include <glas2/backend/current_backend.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/bsp/matrix/concept/bsp_internal_dense_matrix.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>

namespace glas2 {

  template <typename To, typename From>
  typename std::enable_if< is<DenseMatrix,To>::value && is<bsp::BSPInternalDenseMatrix,From>::value, To >::type assign( To to, From const& from ) {
    return assign( current_backend(), to, from ) ;
  }

} // namespace glas2

#endif
