#ifndef glas2_algorithm_norm_1_hpp
#define glas2_algorithm_norm_1_hpp

#include <glas2/concept/external.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>

namespace glas2 {

  template <typename X>
  typename std::enable_if< is<External,X>::value
                         , auto
                         >::type norm_1( X const& x ) -> decltype( norm_1( glas2::wrap(x) ) ) {
    return glas2::norm_1( glas2::wrap(x) ) ;
  }

} // namespace glas2

#endif
