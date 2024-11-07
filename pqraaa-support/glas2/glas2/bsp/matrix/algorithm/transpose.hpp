//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_bsp_matrix_algorithm_transpose_hpp
#define glas2_bsp_matrix_algorithm_transpose_hpp

#include <glas2/matrix/algorithm/transpose.hpp>
#include <glas2/bsp/matrix/type/row_distributed_matrix.hpp>
#include <glas2/bsp/matrix/type/column_distributed_matrix.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>

namespace glas2 {

  template <typename M>
  typename std::enable_if< is<bsp::BSPMatrix,M>::value, typename matrix_transform< transpose_tag, M >::result_type >::type transpose( M m ) {
    return matrix_transform< transpose_tag, M>::apply(m) ;
  }

  template <typename M>
  struct matrix_transform< transpose_tag, M, typename std::enable_if< is<bsp::BSPRowMatrix,M>::value >::type > {
    typedef typename matrix_transform< transpose_tag, typename M::local_type >::result_type local_type ;
    typedef bsp::column_distributed_matrix< local_type, typename M::distribution_type >     result_type ;

    static result_type apply( M m ) {
      return result_type( transpose( m.local() ), m.distribution() ) ;
    }
  } ;

  template <typename M>
  struct matrix_transform< transpose_tag, M, typename std::enable_if< is<bsp::BSPColumnMatrix,M>::value >::type > {
    typedef typename matrix_transform< transpose_tag, typename M::local_type >::result_type local_type ;
    typedef bsp::row_distributed_matrix< local_type, typename M::distribution_type >        result_type ;

    static result_type apply( M m ) {
      return result_type( transpose( m.local() ), m.distribution() ) ;
    }
  } ;

} // namespace glas2

#endif
