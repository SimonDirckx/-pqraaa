#ifndef glas2_algorithm_imag_hpp
#define glas2_algorithm_imag_hpp

#include <glas2/concept/imag_part.hpp>
#include <glas2/type/transformation.hpp>

namespace glas2 {

  template <typename X>
  typename std::enable_if< !std::is_arithmetic<X>::value, transformation< X, imag_part<typename X::value_type> > >::type imag( X x ) {
    return transformation< X, imag_part<typename X::value_type> >( x ) ;
  }

} // namespace glas2

#endif
