#ifndef glas3_bindings_array_dense_array_bind_matrix_hpp
#define glas3_bindings_array_dense_array_bind_matrix_hpp

#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/concept/is.hpp>

#include <glas3/array/dense_array/type/glas2_dense_matrix.hpp>

#include <type_traits>

namespace glas3 {

template <typename M>
typename std::enable_if< glas2::is< glas2::DenseMatrix, M >::value, glas2_dense_matrix< M > >::type
bind_matrix( M& m ) {
	return glas2_dense_matrix< M >( m ) ;
}

} // namespace glas3

#endif
