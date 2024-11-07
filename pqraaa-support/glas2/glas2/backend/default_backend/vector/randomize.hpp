#ifndef glas2_backend_default_backend_vector_randomize_hpp
#define glas2_backend_default_backend_vector_randomize_hpp

#include <glas2/backend/default_backend/default_backend.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/concept/is.hpp>

namespace glas2 {

  template <typename V, typename S>
  typename std::enable_if< is<DenseVector, V>::value, V>::type randomize( default_backend, V& v, S& seed ) {
    for (typename std::decay<V>::type::size_type i=0; i<v.size(); ++i) {
      seed.var_gen( v(i) ) ;
    }
    return v ;
  }

} // namespace glas2

#endif
