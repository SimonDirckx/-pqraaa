//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_backend_default_backend_matrix_algorithm_fill_hpp
#define glas2_backend_default_backend_matrix_algorithm_fill_hpp

#include <glas2/backend/default_backend/default_backend.hpp>
#include <glas2/matrix/algorithm/transpose.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>

namespace glas2 {
/*
  template <typename V, typename T>
  typename std::enable_if< is<DenseMatrix,V>::value, V >::type fill( V v, T const& value ) {
    auto it = v.iterate() ;
    for (typename V::size_type i=0; i<it.size_1(); ++i) {
      for (typename V::size_type j=0; j<it.size_2(); ++j) {
        v[ *it ] = value ;
        it.index_2_pp() ;
      }
      it.index_1_pp() ;
    }
    return v ;
  }
  */

  template <typename M, typename T>
  typename std::enable_if< is<DenseMatrix,M>::value, M >::type fill( current_backend, M m, T const& value ) {
    for (typename M::size_type i=0; i<m.num_rows(); ++i) {
      for (typename M::size_type j=0; j<m.num_columns(); ++j) {
        m(i,j) = value ;
      }
    }
    return m ;
  }

} // namespace glas2

#endif
