#ifndef glas2_matrix_algorithm_outer_prod_hpp
#define glas2_matrix_algorithm_outer_prod_hpp

#include <glas2/concept/is.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/matrix/expression/outer_prod_expression.hpp>
#include <type_traits>

namespace glas2 {

  template <typename V1, typename V2>
  typename std::enable_if< is<DenseVector,V1>::value && is<DenseVector,V2>::value
                         , outer_prod_expression<V1,V2>
                         >::type outer_prod( V1 const& v1, V2 const& v2 ) {
    return outer_prod_expression<V1,V2>( v1, v2 ) ;
  }

} // namespace glas2

#endif
