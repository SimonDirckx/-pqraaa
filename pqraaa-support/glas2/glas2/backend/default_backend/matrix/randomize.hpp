#ifndef glas2_backend_default_backend_matrix_randomize_hpp
#define glas2_backend_default_backend_matrix_randomize_hpp

#include <glas2/backend/default_backend/default_backend.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/concept/is.hpp>
#include <glas2/type/seed.hpp>
#include <type_traits>

namespace glas2 {
/*
  template <typename M, typename S>
  typename std::enable_if< glas2::is<DenseMatrix, M>::value, M>::type randomize( M m, S& seed ) {
    auto it = m.iterate() ;
    for (typename M::size_type i=0; i<it.size_1(); ++i) {
      for (typename M::size_type j=0; j<it.size_2(); ++j) {
        seed.var_gen( m[ (*it) ] ) ;
        it.index_2_pp() ;
      }
      it.index_1_pp() ;
    }
    return m ;
  }
*/

  template <typename M, typename S>
  typename std::enable_if< glas2::is<DenseMatrix, M>::value, M>::type randomize( default_backend, M&& m, S& seed ) {
    typedef typename std::decay<M>::type::size_type size_type ;
    for (size_type i=0; i<m.num_rows(); ++i) {
      for (size_type j=0; j<m.num_columns(); ++j) {
        seed.var_gen( m(i,j) ) ;
      }
    }
    return m ;
  }

} // namespace glas2

#endif
