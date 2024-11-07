#ifndef glas2_algorithm_ops_assign_hpp
#define glas2_algorithm_ops_assign_hpp

#include <glas2/type/pass_reference.hpp>
#include <glas2/backend/current_backend.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>

namespace glas2 { //namespace ops {

  template <typename M, typename E>
  typename std::enable_if< is<Expression,M>::value && pass_reference<M>::pass_by_value, M >::type operator+=( M m, E const& e ) {
    return plus_assign( current_backend(), m, e ) ;
  }

  template <typename M, typename E>
  typename std::enable_if< is<Expression,M>::value && !pass_reference<M>::pass_by_value, typename pass_reference<M>::type >::type operator+=( M& m, E const& e ) {
    typename pass_reference<M>::type m_v = pass_reference<M>()(m) ;
    plus_assign( current_backend(), m_v, e ) ;
    return m_v ;
  }

  template <typename M, typename E>
  typename std::enable_if< is<Expression,M>::value && pass_reference<M>::pass_by_value, M >::type operator-=( M m, E const& e ) {
    return minus_assign( current_backend(), m, e ) ;
  }

  template <typename M, typename E>
  typename std::enable_if< is<Expression,M>::value && !pass_reference<M>::pass_by_value, typename pass_reference<M>::type >::type operator-=( M& m, E const& e ) {
    auto m_v = pass_reference<M>()(m) ;
    return minus_assign( current_backend(), m_v, e ) ;
    return m_v ;
  }

/*
  template <typename M, typename E>
  typename std::enable_if< is<Expression,M>::value, M >::type operator*=( M m, E const& e ) {
    return multiplies_assign( current_backend(), m, e ) ;
  }

  template <typename M, typename E>
  typename std::enable_if< is<Expression,M>::value, M >::type operator/=( M m, E const& e ) {
    return divides_assign( current_backend(), m, e ) ;
  }*/

}// } // namespace glas2::ops

#endif
