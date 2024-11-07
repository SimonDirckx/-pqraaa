#ifndef glas2_backend_openmp_backend_matrix_ops_assign_hpp
#define glas2_backend_openmp_backend_matrix_ops_assign_hpp

#include <glas2/backend/openmp_backend/openmp.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace glas2 {

  template <typename M, typename E>
  typename std::enable_if< is<DenseMatrix,M>::value && is<DenseMatrix,E>::value, M& >::type plus_assign( openmp_backend, M& m, E const& e ) {
    assert( m.num_rows() == e.num_rows() ) ;
    assert( m.num_columns() == e.num_columns() ) ;
#pragma omp parallel for collapse(2)
    for (typename M::size_type i=0; i<m.num_rows(); ++i) {
      for (typename M::size_type j=0; j<m.num_columns(); ++j) {
        m(i,j) += e(i,j) ;
      }
    }
    return m ;
  }

  template <typename M, typename E>
  typename std::enable_if< is<DenseMatrix,M>::value && is<DenseMatrix,E>::value, M& >::type minus_assign( openmp_backend, M& m, E const& e ) {
    assert( m.num_rows() == e.num_rows() ) ;
    assert( m.num_columns() == e.num_columns() ) ;
#pragma omp parallel for collapse(2)
    for (typename M::size_type i=0; i<m.num_rows(); ++i) {
      for (typename M::size_type j=0; j<m.num_columns(); ++j) {
        m(i,j) -= e(i,j) ;
      }
    }
    return m ;
  }

  template <typename M, typename E>
  typename std::enable_if< is<DenseMatrix,M>::value && is<DenseMatrix,E>::value, M& >::type multiplies_assign( openmp_backend, M& m, E const& e ) {
    assert( m.num_rows() == e.num_rows() ) ;
    assert( m.num_columns() == e.num_columns() ) ;
#pragma omp parallel for collapse(2)
    for (typename M::size_type i=0; i<m.num_rows(); ++i) {
      for (typename M::size_type j=0; j<m.num_columns(); ++j) {
        m(i,j) *= e(i,j) ;
      }
    }
    return m ;
  }

  template <typename M, typename E>
  typename std::enable_if< is<DenseMatrix,M>::value && is<Scalar,E>::value, M& >::type multiplies_assign( openmp_backend, M& m, E const& e ) {
#pragma omp parallel for collapse(2)
    for (typename M::size_type i=0; i<m.num_rows(); ++i) {
      for (typename M::size_type j=0; j<m.num_columns(); ++j) {
        m(i,j) *= e ;
      }
    }
    return m ;
  }

  template <typename M, typename E>
  typename std::enable_if< is<DenseMatrix,M>::value && is<Scalar,E>::value, M& >::type divides_assign( openmp_backend, M& m, E const& e ) {
#pragma omp parallel for collapse(2)
    for (typename M::size_type i=0; i<m.num_rows(); ++i) {
      for (typename M::size_type j=0; j<m.num_columns(); ++j) {
        m(i,j) /= e ;
      }
    }
    return m ;
  }

  template <typename M, typename E>
  typename std::enable_if< is<DenseMatrix,M>::value && is<DenseMatrix,E>::value, M& >::type divides_assign( openmp_backend, M& m, E const& e ) {
    assert( m.num_rows() == e.num_rows() ) ;
    assert( m.num_columns() == e.num_columns() ) ;
#pragma omp parallel for collapse(2)
    for (typename M::size_type i=0; i<m.num_rows(); ++i) {
      for (typename M::size_type j=0; j<m.num_columns(); ++j) {
        m(i,j) /= e(i,j) ;
      }
    }
    return m ;
  }

} // namespace glas2

#endif
