#ifndef glas2_algorithm_abs_hpp
#define glas2_algorithm_abs_hpp

#include <glas2/concept/absolute.hpp>
#include <glas2/expression/unary_operation.hpp>
#include <type_traits>

namespace glas2 {

  template <typename X>
  unary_operation< typename std::decay<X>::type, absolute > abs( X&& x ) {
    return unary_operation< typename std::decay<X>::type, absolute >( x ) ;
  }

} // namespace glas2

#endif
