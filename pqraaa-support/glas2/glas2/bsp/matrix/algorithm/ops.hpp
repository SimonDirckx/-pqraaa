#ifndef glas2_bsp_matrix_algorithm_ops_hpp
#define glas2_bsp_matrix_algorithm_ops_hpp

#include <glas2/algorithm/ops.hpp>
#include <glas2/bsp/matrix/expression/binary_operation.hpp>
#include <glas2/bsp/matrix/expression/unary_operation.hpp>

namespace glas2 { namespace bsp {
  using glas2::ops::operator- ;
  using glas2::ops::operator+ ;
  using glas2::ops::operator/ ;
  using glas2::ops::operator* ;
} } // namespace glas2::bsp

#endif
