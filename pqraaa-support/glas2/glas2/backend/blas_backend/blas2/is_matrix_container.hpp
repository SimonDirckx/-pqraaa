#ifndef glas2_backend_blas_backend_matrix_is_matrix_container_hpp
#define glas2_backend_blas_backend_matrix_is_matrix_container_hpp

#include <glas2/matrix/concept/strided_dense_matrix.hpp>
#include <glas2/matrix/concept/contiguous_dense_matrix.hpp>
#include <glas2/matrix/expression/unary_operation.hpp>
#include <glas2/concept/is.hpp>
#include <glas2/concept/conjugate.hpp>

namespace glas2 {

  template <typename M>
  struct is_matrix_container {
    static bool const value = is<StridedDenseMatrix,M>::value || is<ContiguousDenseMatrix,M>::value ;
  } ;

  template <typename M>
  struct is_matrix_container< unary_operation< M, glas2::conjugate > >
  : is_matrix_container< M >
  {} ;

} // namespace glas2

#endif
