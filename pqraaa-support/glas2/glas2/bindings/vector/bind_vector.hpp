#ifndef glas2_bindings_vector_bind_vector_hpp
#define glas2_bindings_vector_bind_vector_hpp

#include <glas2/vector/type/strided_vector.hpp>
#include <boost/numeric/bindings/size.hpp>
#include <boost/numeric/bindings/stride.hpp>
#include <boost/numeric/bindings/begin.hpp>
#include <boost/numeric/bindings/value_type.hpp>
#include <type_traits>

namespace glas2 {

  template <typename V, typename Id, typename EnableIf=void>
  struct bind_vector_selector {
    typedef typename boost::numeric::bindings::value_type<V>::type        value_type ;
    typedef typename boost::numeric::bindings::result_of::size< V >::type size_type ;
    typedef glas2::strided_vector< value_type, size_type >                type ;

    static type apply( Id& v ) {
      return type( boost::numeric::bindings::begin_value(v), boost::numeric::bindings::size(v), boost::numeric::bindings::stride(v) ) ;
    }
  } ;

  template <typename V>
  typename bind_vector_selector< V, V >::type bind_vector( V& v ) {
    return bind_vector_selector< V, V >::apply( v ) ;
  }

/*  template <typename V>
  typename bind_vector_selector< V, const V >::type bind_const_vector( V const& v ) {
    return bind_vector_selector< V, const V >::apply( v ) ;
  }*/

} // namespace glas2

#endif
