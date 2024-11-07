//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_backend_blas_backend_blas2_ops_assign_hpp
#define glas2_backend_blas_backend_blas2_ops_assign_hpp

#include <glas2/backend/blas_backend/blas_backend.hpp>
#include <glas2/backend/blas_backend/blas1/is_vector_container.hpp>
#include <glas2/backend/blas_backend/blas2/is_matrix_container.hpp>
#include <glas2/matrix/expression/matrix_vector_multiply.hpp>
#include <glas2/matrix/expression/outer_prod_expression.hpp>
#include <glas2/vector/expression/unary_operation.hpp>
#include <glas2/bindings/vector/adaptor.hpp>
#include <glas2/bindings/matrix/adaptor.hpp>
#include <boost/numeric/bindings/blas/level2/gemv.hpp>
#include <boost/numeric/bindings/blas/level2/ger.hpp>
#include <boost/numeric/bindings/blas/level2/geru.hpp>
#include <boost/numeric/bindings/blas/level2/gerc.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename W, typename M, typename V>
  typename std::enable_if< is_vector_container<W>::value && is_vector_container<V>::value && is_matrix_container<M>::value, W >::type plus_assign( blas_backend, W w, matrix_vector_multiply<M,V> const& e ) {
    assert( w.size() == e.size() ) ;
    boost::numeric::bindings::blas::gemv( typename W::value_type(1.0), e.matrix(), e.vector(), typename W::value_type(1.0), w ) ;
    return w ;
  }

  template <typename W, typename M, typename V>
  typename std::enable_if< is_vector_container<W>::value && is_vector_container<V>::value && is_matrix_container<M>::value, W >::type minus_assign( blas_backend, W w, matrix_vector_multiply<M,V> const& e ) {
    assert( w.size() == e.size() ) ;
    boost::numeric::bindings::blas::gemv( typename W::value_type(-1.0), e.matrix(), e.vector(), typename W::value_type(1.0), w ) ;
    return w ;
  }

  template <typename W, typename M, typename V>
  typename std::enable_if< is_vector_container<W>::value && is_vector_container<V>::value && is_matrix_container<M>::value && std::is_arithmetic<typename M::value_type>::value, M& >::type plus_assign( blas_backend, M& m, outer_prod_expression<V,W> const& e ) {
    assert( m.num_rows() == e.num_rows() ) ;
    assert( m.num_columns() == e.num_columns() ) ;
    boost::numeric::bindings::blas::ger( typename M::value_type(1.0), e.vector1(), e.vector2(), m ) ;
    return m ;
  }

  template <typename W, typename M, typename V>
  typename std::enable_if< is_vector_container<W>::value && is_vector_container<V>::value && is_matrix_container<M>::value && std::is_arithmetic<typename M::value_type>::value, M& >::type minus_assign( blas_backend, M& m, outer_prod_expression<V,W> const& e ) {
    assert( m.num_rows() == e.num_rows() ) ;
    assert( m.num_columns() == e.num_columns() ) ;
    boost::numeric::bindings::blas::ger( typename M::value_type(-1.0), e.vector1(), e.vector2(), m ) ;
    return m ;
  }

  template <typename W, typename M, typename V>
  typename std::enable_if< is_vector_container<W>::value && is_vector_container<V>::value && is_matrix_container<M>::value && !std::is_arithmetic<typename M::value_type>::value, M& >::type plus_assign( blas_backend, M& m, outer_prod_expression<V,W> const& e ) {
    assert( m.num_rows() == e.num_rows() ) ;
    assert( m.num_columns() == e.num_columns() ) ;
    boost::numeric::bindings::blas::geru( typename M::value_type(1.0), e.vector1(), e.vector2(), m ) ;
    return m ;
  }

  template <typename W, typename M, typename V>
  typename std::enable_if< is_vector_container<W>::value && is_vector_container<V>::value && is_matrix_container<M>::value && !std::is_arithmetic<typename M::value_type>::value, M& >::type minus_assign( blas_backend, M& m, outer_prod_expression<V,W> const& e ) {
    assert( m.num_rows() == e.num_rows() ) ;
    assert( m.num_columns() == e.num_columns() ) ;
    boost::numeric::bindings::blas::geru( typename M::value_type(-1.0), e.vector1(), e.vector2(), m ) ;
    return m ;
  }

  template <typename W, typename M, typename V>
  typename std::enable_if< is_vector_container<W>::value && is_vector_container<V>::value && is_matrix_container<M>::value && !std::is_arithmetic<typename M::value_type>::value, M& >::type minus_assign( blas_backend, M& m, outer_prod_expression<V, unary_operation<W,conjugate> > const& e ) {
    assert( m.num_rows() == e.num_rows() ) ;
    assert( m.num_columns() == e.num_columns() ) ;
    boost::numeric::bindings::blas::gerc( typename M::value_type(-1.0), e.vector1(), e.vector2().vector(), m ) ;
    return m ;
  }

} // namespace glas2

#endif
