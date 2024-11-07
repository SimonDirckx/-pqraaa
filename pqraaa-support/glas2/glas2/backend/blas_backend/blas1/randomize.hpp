#ifndef glas2_backend_blas_backend_blas1_randomize_hpp
#define glas2_backend_blas_backend_blas1_randomize_hpp

#include <glas2/backend/blas_backend/blas_backend.hpp>
#include <glas2/backend/default_backend/vector/randomize.hpp>

namespace glas2 {

  template <typename V, typename S>
  typename std::enable_if< is<DenseVector, V>::value, V>::type randomize( blas_backend, V v, S& seed ) {
    return randomize( default_backend(), v, seed ) ;
  }

} // namespace glas2

#endif
