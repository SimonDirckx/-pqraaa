//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_sparse_algorithm_upper_hpp
#define glas2_sparse_algorithm_upper_hpp

#include <glas2/sparse/concept/sparse_matrix.hpp>
#include <glas2/sparse/view/sparse_upper_view.hpp>
#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <type_traits>
#include <cmath>

namespace glas2 {

  template <typename M>
  typename std::enable_if< is< SparseMatrix, M >::value, sparse_upper_view<M> >::type upper( M m ) {
    return sparse_upper_view<M>( m ) ;
  }

} // namespace glas2

#endif
