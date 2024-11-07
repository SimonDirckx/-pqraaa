#ifndef glas2_expression_binary_operation_hpp
#define glas2_expression_binary_operation_hpp

#include <glas2/concept/is.hpp>
#include <type_traits>

namespace glas2 {

  template <typename E1, typename E2, typename Op, typename EnableIf=void>
  struct binary_operation {
  } ;

} // namespace glas2

#endif
