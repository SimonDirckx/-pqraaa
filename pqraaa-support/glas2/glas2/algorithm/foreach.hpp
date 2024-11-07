#ifndef glas2_algorithm_foreach_hpp
#define glas2_algorithm_foreach_hpp

#include <glas2/expression/binary_operation.hpp>
#include <glas2/expression/unary_operation.hpp>
#include <glas2/concept/is.hpp>
#include <glas2/concept/expression.hpp>
#include <type_traits>

namespace glas2 {

  template <typename A, typename F>
  unary_operation<A,F> foreach( A const& a, F const& f ) {
    return unary_operation<A,F>( a, f ) ;
  }

} // namespace glas2

#endif
