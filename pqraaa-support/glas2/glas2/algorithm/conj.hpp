#ifndef glas2_algorithm_conj_hpp
#define glas2_algorithm_conj_hpp

#include <glas2/concept/conjugate.hpp>
#include <glas2/expression/unary_operation.hpp>

namespace glas2 {

  template <typename X>
  typename std::enable_if< !std::is_arithmetic<X>::value, unary_operation< X, conjugate > >::type conj( X const& x ) {
    return unary_operation< X, conjugate >( x ) ;
  }

} // namespace glas2

#endif
