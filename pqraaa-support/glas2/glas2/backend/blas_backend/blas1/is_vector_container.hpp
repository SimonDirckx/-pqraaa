#ifndef glas2_backend_blas_backend_vector_is_vector_container_hpp
#define glas2_backend_blas_backend_vector_is_vector_container_hpp

#include <glas2/vector/concept/strided_dense_vector.hpp>
#include <glas2/vector/concept/contiguous_dense_vector.hpp>
#include <glas2/concept/is.hpp>

namespace glas2 {

  template <typename V>
  struct is_vector_container {
    static bool const value = is<StridedDenseVector,V>::value || is<ContiguousDenseVector,V>::value ;
  } ;

} // namespace glas2

#endif
