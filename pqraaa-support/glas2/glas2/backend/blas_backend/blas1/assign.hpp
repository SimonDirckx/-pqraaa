//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_backend_blas_backend_blas1_assign_hpp
#define glas2_backend_blas_backend_blas1_assign_hpp

#include <glas2/backend/blas_backend/blas_backend.hpp>
#include <glas2/backend/blas_backend/blas1/is_vector_container.hpp>
#include <glas2/vector/concept/strided_dense_vector.hpp>
#include <glas2/vector/concept/contiguous_dense_vector.hpp>
#include <glas2/vector/expression/binary_operation.hpp>
#include <glas2/vector/expression/unary_operation.hpp>
#include <glas2/concept/negate.hpp>
#include <glas2/bindings/vector/adaptor.hpp>
#include <boost/numeric/bindings/blas/level1/copy.hpp>
#include <boost/numeric/bindings/blas/level1/axpy.hpp>
#include <boost/numeric/bindings/blas/level1/scal.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename To, typename From>
  typename std::enable_if< is_vector_container<To>::value && is_vector_container<From>::value, To >::type assign( blas_backend, To to, From const& from ) {
    boost::numeric::bindings::blas::copy( from, to ) ;
    return to ;
  }

  template <typename V, typename E, typename S>
  typename std::enable_if< is_vector_container<V>::value && is_vector_container<E>::value && is<Scalar,S>::value, V >::type assign( blas_backend, V v, binary_operation<S,E,multiplies> const& e ) {
    assert( v.size() == e.size() ) ;
    boost::numeric::bindings::blas::scal( v, typename V::value_type(0.0) ) ;
    boost::numeric::bindings::blas::axpy( e.scalar(), e.vector(), v ) ;
    return v ;
  }

  template <typename V, typename E, typename S>
  typename std::enable_if< is_vector_container<V>::value && is_vector_container<E>::value && is<Scalar,S>::value, V >::type assign( blas_backend, V v, binary_operation<E,S,multiplies> const& e ) {
    assert( v.size() == e.size() ) ;
    boost::numeric::bindings::blas::scal( v, typename V::value_type(0.0) ) ;
    boost::numeric::bindings::blas::axpy( e.scalar(), e.vector(), v ) ;
    return v ;
  }

  template <typename V, typename E>
  typename std::enable_if< is_vector_container<V>::value && is_vector_container<E>::value, V >::type assign( blas_backend, V v, unary_operation<E,negate> const& e ) {
    assert( v.size() == e.size() ) ;
    boost::numeric::bindings::blas::scal( v, typename V::value_type(0.0) ) ;
    boost::numeric::bindings::blas::axpy( typename E::value_type(-1.0), e.vector(), v ) ;
    return v ;
  }

} // namespace glas2

#endif
