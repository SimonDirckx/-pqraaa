//  (C) Copyright Karl Meerbergen 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_sparse_algorithm_norm_fro_hpp
#define glas2_sparse_algorithm_norm_fro_hpp

#include <glas2/vector/algorithm/norm_2.hpp>
#include <glas2/sparse/concept/compressed_sparse_matrix.hpp>
#include <glas2/sparse/concept/coordinate_sparse_matrix.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cmath>

namespace glas2 {

  template <typename M>
  typename std::enable_if< ::glas2::is< ::glas2::CompressedSparseMatrix, M >::value, decltype( std::abs(typename M::value_type() ) ) >::type norm_fro( M const& m ) {
    return norm_2( m.data() ) ;
  }

  template <typename M>
  typename std::enable_if< ::glas2::is< ::glas2::CoordinateSparseMatrix, M >::value, decltype( std::abs(typename M::value_type() ) ) >::type norm_fro( M const& m ) {
    return norm_2( m.data() ) ;
  }
}

#endif
