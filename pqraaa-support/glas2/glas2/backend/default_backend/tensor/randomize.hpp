#ifndef glas2_backend_default_backend_tensor_randomize_hpp
#define glas2_backend_default_backend_tensor_randomize_hpp

#include <glas2/backend/default_backend/default_backend.hpp>
#include <glas2/tensor/concept/dense_tensor.hpp>
#include <glas2/concept/is.hpp>
#include <glas2/type/seed.hpp>
#include <algorithm>
#include <type_traits>
#include <vector>

namespace glas2 {

  template <typename X, typename S>
  typename std::enable_if< glas2::is<DenseTensor, X>::value, X>::type randomize( default_backend, X x, S& seed ) {
    typedef typename X::size_type size_type ;
    std::vector<size_type> index( x.order()+1 ) ;
    std::fill( index.begin(), index.end(), 0 ) ;

    size_type nnz = prod( x.shape() ) ;
    for (size_type i=0; i<nnz; ++i) {
      seed.var_gen( x[i] ) ;
      index[0]++ ; int k=0; while ( k<x.order() && index[k]==x.shape()(k) ) { index[k] = 0; ++k; index[k]++ ; }
    }
    return x ;
  }

} // namespace glas2

#endif
