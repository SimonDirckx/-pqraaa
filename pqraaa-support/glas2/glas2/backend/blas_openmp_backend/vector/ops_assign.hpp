#ifndef glas2_backend_blas_openmp_vector_algorithm_ops_assign_hpp
#define glas2_backend_blas_openmp_vector_algorithm_ops_assign_hpp

#include <glas2/backend/blas_backend/blas1.hpp>
#include <glas2/backend/blas_backend/blas2.hpp>
#include <glas2/backend/blas_openmp_backend/blas_openmp.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cassert>

namespace glas2 {

  template <typename V, typename E>
  typename std::enable_if< is<DenseVector,V>::value && is<DenseVector,E>::value, V& >::type plus_assign( blas_openmp, V& v, E const& e ) {
    assert( v.size() == e.size() ) ;
#ifdef _OPENMP
    int n_threads = omp_get_num_threads() ;
#else
    int n_threads = 1 ;
#endif
    typename V::size_type step = (v.size()+n_threads-1) / n_threads ;
#pragma omp parallel for
    for (int i=0; i<n_threads; ++i) {
      glas2::range r( step*i, std::min(v.size(),step*i+step) ) ;
      auto to_r = v( r ) ;
      plus_assign( blas_backend(), to_r, e(r) ) ;
    }
    return v ;
  }

  template <typename V, typename E>
  typename std::enable_if< is<DenseVector,V>::value && is<DenseVector,E>::value, V& >::type minus_assign( blas_openmp, V& v, E const& e ) {
    assert( v.size() == e.size() ) ;
#ifdef _OPENMP
    int n_threads = omp_get_num_threads() ;
#else
    int n_threads = 1 ;
#endif
    typename V::size_type step = (v.size()+n_threads-1) / n_threads ;
#pragma omp parallel for
    for (int i=0; i<n_threads; ++i) {
      glas2::range r( step*i, std::min(v.size(),step*i+step) ) ;
      auto to_r = v( r ) ;
      minus_assign( blas_backend(), to_r, e(r) ) ;
    }
    return v ;
  }

  template <typename V, typename E>
  typename std::enable_if< is<DenseVector,V>::value && is<Scalar,E>::value, V& >::type multiplies_assign( blas_openmp, V& v, E const& e ) {
#ifdef _OPENMP
    int n_threads = omp_get_num_threads() ;
#else
    int n_threads = 1 ;
#endif
    typename V::size_type step = (v.size()+n_threads-1) / n_threads ;
#pragma omp parallel for
    for (int i=0; i<n_threads; ++i) {
      glas2::range r( step*i, std::min(v.size(),step*i+step) ) ;
      auto to_r = v( r ) ;
      multiplies_assign( blas_backend(), to_r, e(r) ) ;
    }
    return v ;
  }

} // namespace glas2

#endif
