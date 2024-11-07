#ifndef glas2_backend_openmp_vector_algorithm_ops_assign_hpp
#define glas2_backend_openmp_vector_algorithm_ops_assign_hpp

#include <glas2/backend/openmp_backend/openmp.hpp>
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
  typename std::enable_if< is<DenseVector,V>::value && is<DenseVector,E>::value, V& >::type plus_assign( openmp, V& v, E const& e ) {
    assert( v.size() == e.size() ) ;
#pragma omp parallel for
    for (typename V::size_type i=0; i<v.size(); ++i) {
      v(i) += e(i) ;
    }
    return v ;
  }

  template <typename V, typename E>
  typename std::enable_if< is<DenseVector,V>::value && is<Scalar,E>::value, V& >::type plus_assign( openmp, V& v, E const& e ) {
#pragma omp parallel for
    for (typename V::size_type i=0; i<v.size(); ++i) {
      v(i) += e ;
    }
    return v ;
  }

  template <typename V, typename E>
  typename std::enable_if< is<DenseVector,V>::value && is<DenseVector,E>::value, V& >::type minus_assign( openmp, V& v, E const& e ) {
    assert( v.size() == e.size() ) ;
#pragma omp parallel for
    for (typename V::size_type i=0; i<v.size(); ++i) {
      v(i) -= e(i) ;
    }
    return v ;
  }

  template <typename V, typename E>
  typename std::enable_if< is<DenseVector,V>::value && is<Scalar,E>::value, V& >::type minus_assign( openmp, V& v, E const& e ) {
#pragma omp parallel for
    for (typename V::size_type i=0; i<v.size(); ++i) {
      v(i) -= e ;
    }
    return v ;
  }

  template <typename V, typename E>
  typename std::enable_if< is<DenseVector,V>::value && is<Scalar,E>::value, V& >::type multiplies_assign( openmp, V& v, E const& e ) {
#pragma omp parallel for
    for (typename V::size_type i=0; i<v.size(); ++i) {
      v(i) *= e ;
    }
    return v ;
  }

  template <typename V, typename E>
  typename std::enable_if< is<DenseVector,V>::value && is<Scalar,E>::value, V& >::type divides_assign( openmp, V& v, E const& e ) {
#pragma omp parallel for
    for (typename V::size_type i=0; i<v.size(); ++i) {
      v(i) /= e ;
    }
    return v ;
  }

} // namespace glas2

#endif
