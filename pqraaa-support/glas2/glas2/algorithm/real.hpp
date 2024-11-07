#ifndef glas2_algorithm_real_hpp
#define glas2_algorithm_real_hpp

#include <glas2/concept/real_part.hpp>
#include <glas2/type/transformation.hpp>

namespace glas2 {

  template <typename X>
  typename std::enable_if< !std::is_arithmetic<X>::value, transformation< X, real_part > >::type real( X x ) {
    return transformation< X, real_part >( x ) ;
  }

} // namespace glas2

#endif
