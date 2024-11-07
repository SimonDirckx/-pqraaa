#ifndef glas2_backend_default_backend_vector_algorithm_ops_assign_hpp
#define glas2_backend_default_backend_vector_algorithm_ops_assign_hpp

#include <glas2/backend/default_backend/default_backend.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename V, typename E>
  typename std::enable_if< is<DenseVector,V>::value && is<DenseVector,E>::value, V& >::type plus_assign( default_backend, V& v, E const& e ) {
    assert( v.size() == e.size() ) ;
    for (typename V::size_type i=0; i<v.size(); ++i) {
      v(i) += e(i) ;
    }
    return v ;
  }

  template <typename V, typename E>
  typename std::enable_if< is<DenseVector,V>::value && is<Scalar,E>::value, V& >::type plus_assign( default_backend, V& v, E const& e ) {
    for (typename V::size_type i=0; i<v.size(); ++i) {
      v(i) += e ;
    }
    return v ;
  }

  template <typename V, typename E>
  typename std::enable_if< is<DenseVector,V>::value && is<DenseVector,E>::value, V& >::type minus_assign( default_backend, V& v, E const& e ) {
    assert( v.size() == e.size() ) ;
    for (typename V::size_type i=0; i<v.size(); ++i) {
      v(i) -= e(i) ;
    }
    return v ;
  }

  template <typename V, typename E>
  typename std::enable_if< is<DenseVector,V>::value && is<Scalar,E>::value, V& >::type minus_assign( default_backend, V& v, E const& e ) {
    for (typename V::size_type i=0; i<v.size(); ++i) {
      v(i) -= e ;
    }
    return v ;
  }

  template <typename V, typename E>
  typename std::enable_if< is<DenseVector,V>::value && is<Scalar,E>::value, V& >::type multiplies_assign( default_backend, V& v, E const& e ) {
    for (typename V::size_type i=0; i<v.size(); ++i) {
      v(i) *= e ;
    }
    return v ;
  }

  template <typename V, typename E>
  typename std::enable_if< is<DenseVector,V>::value && is<Scalar,E>::value, V& >::type divides_assign( default_backend, V& v, E const& e ) {
    for (typename V::size_type i=0; i<v.size(); ++i) {
      v(i) /= e ;
    }
    return v ;
  }

  template <typename V, typename E>
  typename std::enable_if< is<DenseVector,V>::value && is<DenseVector,E>::value, V& >::type divides_assign( default_backend, V& v, E const& e ) {
    assert( v.size()==e.size() ) ;
    for (typename V::size_type i=0; i<v.size(); ++i) {
      v(i) /= e(i) ;
    }
    return v ;
  }

} // namespace glas2

#endif
