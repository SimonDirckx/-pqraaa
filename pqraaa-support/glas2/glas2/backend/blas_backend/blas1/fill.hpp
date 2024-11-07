#ifndef glas2_backend_blas_backend_blas1_fill_hpp
#define glas2_backend_blas_backend_blas1_fill_hpp

#include <glas2/backend/blas_backend/blas_backend.hpp>
#include <glas2/backend/default_backend/vector/fill.hpp>

namespace glas2 {

  template <typename V, typename E>
  typename std::enable_if< is<DenseVector, V>::value, V>::type fill( blas_backend, V v, E const& e ) {
    return fill( default_backend(), v, e ) ;
  }

} // namespace glas2

#endif
