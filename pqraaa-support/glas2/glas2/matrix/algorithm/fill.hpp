//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_algorithm_fill_hpp
#define glas2_matrix_algorithm_fill_hpp

#include <glas2/backend/default_backend/matrix/fill.hpp>
#include <glas2/backend/current_backend.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>

namespace glas2 {

  template <typename M, typename T>
  typename std::enable_if< is<DenseMatrix,M>::value, M >::type fill( M m, T const& value ) {
    return fill( current_backend(), m, value ) ;
  }

} // namespace glas2

#endif
