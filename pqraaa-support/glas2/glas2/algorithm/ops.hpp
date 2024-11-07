#ifndef glas2_algorithm_ops_hpp
#define glas2_algorithm_ops_hpp

#include <glas2/expression/binary_operation.hpp>
#include <glas2/expression/unary_operation.hpp>
#include <glas2/concept/divides.hpp>
#include <glas2/concept/is.hpp>
#include <glas2/concept/expression.hpp>
#include <glas2/concept/multiplies.hpp>
#include <glas2/concept/minus.hpp>
#include <glas2/concept/negate.hpp>
#include <glas2/concept/plus.hpp>
#include <glas2/concept/smaller.hpp>
#include <glas2/concept/smaller_equal.hpp>
#include <glas2/concept/equal.hpp>
#include <type_traits>

namespace glas2 { namespace ops {

  template <typename A1, typename A2>
  typename std::enable_if< is<Expression,A1>::value || is<Expression,A2>::value, binary_operation<A1,A2,multiplies> >::type operator*( A1 const& a1, A2 const& a2 ) {
    return binary_operation<A1,A2,multiplies>( a1, a2 ) ;
  }

  template <typename A1, typename A2>
  typename std::enable_if< is<Expression,A1>::value || is<Expression,A2>::value, binary_operation<A1,A2,plus> >::type operator+( A1 const& a1, A2 const& a2 ) {
    return binary_operation<A1,A2,plus>( a1, a2 ) ;
  }

  template <typename A1, typename A2>
  typename std::enable_if< is<Expression,A1>::value || is<Expression,A2>::value, binary_operation<A1,A2,minus> >::type operator-( A1 const& a1, A2 const& a2 ) {
    return binary_operation<A1,A2,minus>( a1, a2 ) ;
  }

  template <typename A1, typename A2>
  typename std::enable_if< is<Expression,A1>::value || is<Expression,A2>::value, binary_operation<A1,A2,divides> >::type operator/( A1 const& a1, A2 const& a2 ) {
    return binary_operation<A1,A2,divides>( a1, a2 ) ;
  }

  template <typename A1, typename A2>
  typename std::enable_if< is<Expression,A1>::value || is<Expression,A2>::value, binary_operation<A1,A2,smaller> >::type operator<( A1 const& a1, A2 const& a2 ) {
    return binary_operation<A1,A2,smaller>( a1, a2 ) ;
  }

  template <typename A1, typename A2>
  typename std::enable_if< is<Expression,A1>::value || is<Expression,A2>::value, binary_operation<A1,A2,smaller_equal> >::type operator<=( A1 const& a1, A2 const& a2 ) {
    return binary_operation<A1,A2,smaller_equal>( a1, a2 ) ;
  }

  template <typename A1, typename A2>
  typename std::enable_if< is<Expression,A1>::value || is<Expression,A2>::value, binary_operation<A1,A2,equal> >::type operator==( A1 const& a1, A2 const& a2 ) {
    return binary_operation<A1,A2,equal>( a1, a2 ) ;
  }

  template <typename A>
  typename std::enable_if< is<Expression,A>::value, unary_operation<A,negate > >::type operator-( A const& a ) {
    return unary_operation<A,negate >( a ) ;
  }

} } // namespace glas2::ops

#endif
