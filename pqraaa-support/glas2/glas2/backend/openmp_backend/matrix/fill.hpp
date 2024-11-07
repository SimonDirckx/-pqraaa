//  (C) Copyright Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_backend_openmp_backend_matrix_algorithm_fill_hpp
#define glas2_backend_openmp_backend_matrix_algorithm_fill_hpp

#include <glas2/backend/openmp_backend/openmp.hpp>
#include <glas2/matrix/algorithm/transpose.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/concept/is.hpp>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cassert>
#include <type_traits>

namespace glas2 {

  template <typename M, typename T>
  typename std::enable_if< is<DenseMatrix,M>::value, M >::type fill( current_backend, M m, T const& value ) {
#pragma omp parallel for collapse(2)
    for (typename M::size_type i=0; i<m.num_rows(); ++i) {
      for (typename M::size_type j=0; j<m.num_columns(); ++j) {
        m(i,j) = value ;
      }
    }
    return m ;
  }

} // namespace glas2

#endif
