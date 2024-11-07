#ifndef glas2_backend_blas_backend_blas1_ops_assign_hpp
#define glas2_backend_blas_backend_blas1_ops_assign_hpp

#include <glas2/backend/blas_backend/blas_backend.hpp>
#include <glas2/backend/blas_backend/blas1/is_vector_container.hpp>
#include <glas2/bindings/vector/adaptor.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/concept/multiplies.hpp>
#include <glas2/concept/divides.hpp>
#include <glas2/vector/concept/strided_dense_vector.hpp>
#include <glas2/vector/concept/contiguous_dense_vector.hpp>
#include <glas2/vector/expression/binary_operation.hpp>
#include <glas2/concept/is.hpp>
#include <boost/numeric/bindings/blas/level1/axpy.hpp>
#include <boost/numeric/bindings/blas/level1/scal.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename V, typename E>
  typename std::enable_if< is_vector_container<V>::value && is_vector_container<E>::value, V >::type plus_assign( blas_backend, V v, E const& e ) {
    assert( v.size() == e.size() ) ;
    boost::numeric::bindings::blas::axpy( typename V::value_type(1.0), e, v ) ;
    return v ;
  }

  template <typename V, typename E>
  typename std::enable_if< is_vector_container<V>::value && is_vector_container<E>::value, V >::type minus_assign( blas_backend, V v, E const& e ) {
    assert( v.size() == e.size() ) ;
    boost::numeric::bindings::blas::axpy( typename V::value_type(-1.0), e, v ) ;
    return v ;
  }

  template <typename V, typename E, typename S>
  typename std::enable_if< is_vector_container<V>::value && is_vector_container<E>::value && is<Scalar,S>::value, V >::type plus_assign( blas_backend, V v, binary_operation<S,E,multiplies> const& e ) {
    assert( v.size() == e.size() ) ;
    boost::numeric::bindings::blas::axpy( e.scalar(), e.vector(), v ) ;
    return v ;
  }

  template <typename V, typename E, typename S>
  typename std::enable_if< is_vector_container<V>::value && is_vector_container<E>::value && is<Scalar,S>::value, V >::type minus_assign( blas_backend, V v, binary_operation<S,E,multiplies> const& e ) {
    assert( v.size() == e.size() ) ;
    boost::numeric::bindings::blas::axpy( -e.scalar(), e.vector(), v ) ;
    return v ;
  }

  template <typename V, typename E, typename S>
  typename std::enable_if< is_vector_container<V>::value && is_vector_container<E>::value && is<Scalar,S>::value, V >::type plus_assign( blas_backend, V v, binary_operation<E,S,multiplies> const& e ) {
    assert( v.size() == e.size() ) ;
    boost::numeric::bindings::blas::axpy( e.scalar(), e.vector(), v ) ;
    return v ;
  }

  template <typename V, typename E, typename S>
  typename std::enable_if< is_vector_container<V>::value && is_vector_container<E>::value && is<Scalar,S>::value, V >::type minus_assign( blas_backend, V v, binary_operation<E,S,multiplies> const& e ) {
    assert( v.size() == e.size() ) ;
    boost::numeric::bindings::blas::axpy( -e.scalar(), e.vector(), v ) ;
    return v ;
  }

  template <typename V, typename E>
  typename std::enable_if< is_vector_container<V>::value && is<Scalar,E>::value, V >::type multiplies_assign( blas_backend, V v, E const& e ) {
    boost::numeric::bindings::blas::scal( e, v ) ;
    return v ;
  }

  template <typename V, typename E>
  typename std::enable_if< is_vector_container<V>::value && is<Scalar,E>::value, V >::type divides_assign( blas_backend, V v, E const& e ) {
    boost::numeric::bindings::blas::scal( 1./e, v ) ;
    return v ;
  }

} // namespace glas2

#endif
