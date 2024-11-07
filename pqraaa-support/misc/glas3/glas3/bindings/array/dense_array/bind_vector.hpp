#ifndef glas3_bindings_array_dense_array_bind_vector_hpp
#define glas3_bindings_array_dense_array_bind_vector_hpp

#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/concept/is.hpp>

#include <glas3/array/dense_array/type/glas2_dense_vector.hpp>

#include <type_traits>

namespace glas3 {

template <typename V>
typename std::enable_if< glas2::is< glas2::DenseVector, V >::value, glas2_dense_vector< V > >::type
bind_vector( V& m ) {
	return glas2_dense_vector< V >( m ) ;
}

} // namespace glas3

#endif
