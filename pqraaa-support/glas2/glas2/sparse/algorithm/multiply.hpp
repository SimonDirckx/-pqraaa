//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_sparse_algorithm_multiply_hpp
#define glas2_sparse_algorithm_multiply_hpp

#include <glas2/concept/is.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/sparse/expression/sparse_matrix_vector_multiply.hpp>
#include <glas2/sparse/expression/sparse_matrix_multiply.hpp>
#include <type_traits>

namespace glas2 {

  template <typename M, typename V>
  typename std::enable_if< is<SparseMatrix,M>::value && is<DenseVector,V>::value
                         , sparse_matrix_vector_multiply<M,V>
                         >::type multiply( M const& m, V const& v ) {
    return sparse_matrix_vector_multiply<M,V>( m, v ) ;
  }

  template <typename M1, typename M2>
  typename std::enable_if< is<SparseMatrix,M1>::value && is<DenseMatrix,M2>::value
                         , sparse_matrix_multiply<M1,M2>
                         >::type multiply( M1 const& m1, M2 const& m2 ) {
    return sparse_matrix_multiply<M1,M2>( m1, m2 ) ;
  }

} // namespace glas2

#endif
