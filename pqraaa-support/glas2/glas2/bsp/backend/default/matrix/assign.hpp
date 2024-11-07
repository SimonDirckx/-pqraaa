//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_bsp_backend_default_matrix_assign_hpp
#define glas2_bsp_backend_default_matrix_assign_hpp

#include <glas2/bsp/backend/default/default_backend.hpp>
#include <glas2/matrix/algorithm/multiply.hpp>
#include <glas2/backend/current_backend.hpp>
#include <glas2/bsp/matrix/concept/bsp_matrix.hpp>
#include <glas2/bsp/matrix/concept/bsp_internal_dense_matrix.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 { namespace bsp {

  template <typename To, typename From>
  typename std::enable_if< ::glas2::is<BSPMatrix,To>::value && ::glas2::is<BSPMatrix,From>::value, To >::type assign( current_backend, To to, From const& from ) {
    // We have to test whether there is no BSP needed in evaluating from
    assert( to.distribution()==from.distribution() ) ;
    to.local() = from.local() ;
    return to ;
  }

} } // namespace glas2::bsp

namespace glas2 {

  template <typename To, typename From>
  typename std::enable_if< ::glas2::is<DenseMatrix,To>::value && ::glas2::is<bsp::BSPInternalDenseMatrix,From>::value, To >::type assign( current_backend, To to, From const& from ) {
    to = multiply( from.matrix1().local(), from.matrix2().local() ) ;
    // Now Allreduce to
    return to ;
  }

} // namespace glas2

#endif
