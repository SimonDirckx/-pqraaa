#ifndef glas2_backend_default_backend_matrix_ops_assign_hpp
#define glas2_backend_default_backend_matrix_ops_assign_hpp

#include <glas2/backend/default_backend/default_backend.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename M, typename E>
  typename std::enable_if< is<DenseMatrix,M>::value && is<DenseMatrix,E>::value, M& >::type plus_assign( default_backend, M& m, E const& e ) {
    assert( m.num_rows() == e.num_rows() ) ;
    assert( m.num_columns() == e.num_columns() ) ;
    for (typename M::size_type i=0; i<m.num_rows(); ++i) {
      for (typename M::size_type j=0; j<m.num_columns(); ++j) {
        m(i,j) += e(i,j) ;
      }
    }
    return m ;
  }

  template <typename M, typename E>
  typename std::enable_if< is<DenseMatrix,M>::value && is<DenseMatrix,E>::value, M& >::type minus_assign( default_backend, M& m, E const& e ) {
    assert( m.num_rows() == e.num_rows() ) ;
    assert( m.num_columns() == e.num_columns() ) ;
    for (typename M::size_type i=0; i<m.num_rows(); ++i) {
      for (typename M::size_type j=0; j<m.num_columns(); ++j) {
        m(i,j) -= e(i,j) ;
      }
    }
    return m ;
  }

  template <typename M, typename E>
  typename std::enable_if< is<DenseMatrix,M>::value && is<DenseMatrix,E>::value, M& >::type multiplies_assign( default_backend, M& m, E const& e ) {
    assert( m.num_rows() == e.num_rows() ) ;
    assert( m.num_columns() == e.num_columns() ) ;
    for (typename M::size_type i=0; i<m.num_rows(); ++i) {
      for (typename M::size_type j=0; j<m.num_columns(); ++j) {
        m(i,j) *= e(i,j) ;
      }
    }
    return m ;
  }

  template <typename M, typename E>
  typename std::enable_if< is<DenseMatrix,M>::value && is<Scalar,E>::value, M& >::type multiplies_assign( default_backend, M& m, E const& e ) {
    for (typename M::size_type i=0; i<m.num_rows(); ++i) {
      for (typename M::size_type j=0; j<m.num_columns(); ++j) {
        m(i,j) *= e ;
      }
    }
    return m ;
  }

  template <typename M, typename E>
  typename std::enable_if< is<DenseMatrix,M>::value && is<Scalar,E>::value, M& >::type divides_assign( default_backend, M& m, E const& e ) {
    for (typename M::size_type i=0; i<m.num_rows(); ++i) {
      for (typename M::size_type j=0; j<m.num_columns(); ++j) {
        m(i,j) /= e ;
      }
    }
    return m ;
  }

  template <typename M, typename E>
  typename std::enable_if< is<DenseMatrix,M>::value && is<DenseMatrix,E>::value, M& >::type divides_assign( default_backend, M& m, E const& e ) {
    assert( m.num_rows() == e.num_rows() ) ;
    assert( m.num_columns() == e.num_columns() ) ;
    for (typename M::size_type i=0; i<m.num_rows(); ++i) {
      for (typename M::size_type j=0; j<m.num_columns(); ++j) {
        m(i,j) /= e(i,j) ;
      }
    }
    return m ;
  }

} // namespace glas2

#endif
