#ifndef glas2_expression_unary_operation_hpp
#define glas2_expression_unary_operation_hpp

#include <glas2/concept/is.hpp>
#include <type_traits>

namespace glas2 {

  template <typename E, typename Op, typename EnableIf=void>
  struct unary_operation {
  } ;

} // namespace glas2

#endif
