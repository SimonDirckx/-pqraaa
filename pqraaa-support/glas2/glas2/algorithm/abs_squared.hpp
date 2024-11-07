#ifndef glas2_algorithm_abs_squared_hpp
#define glas2_algorithm_abs_squared_hpp

#include <glas2/concept/absolute_squared.hpp>
#include <glas2/concept/is_arithmetic.hpp>
#include <glas2/expression/unary_operation.hpp>
#include <type_traits>

namespace glas2 {

  template <typename X>
  typename std::enable_if< !glas2::is_arithmetic<typename std::decay<X>::type>::value, unary_operation< typename std::decay<X>::type, absolute_squared > >::type abs_squared( X&& x ) {
    return unary_operation< typename std::decay<X>::type, absolute_squared >( x ) ;
  }

} // namespace glas2

#endif
