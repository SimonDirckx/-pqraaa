//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_algorithm_vec_hpp
#define glas2_matrix_algorithm_vec_hpp

#include <glas2/matrix/type/contiguous_matrix.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/matrix/concept/matrix_transform.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>

namespace glas2 {

  struct vec_tag {} ;

  template <typename M>
  typename std::enable_if< is<DenseMatrix,M>::value, typename matrix_transform< vec_tag, M >::result_type >::type vec( M m ) {
    //static_assert( std::is_same< column_major, typename M::orientation >::value, "vec: matrix must be column_major" ) ;
    return matrix_transform< vec_tag, M>::apply(m) ;
  }

  template <typename T, typename S>
  struct matrix_transform< vec_tag, contiguous_matrix<T,S,column_major> > {
    typedef contiguous_vector<T,S> result_type ;

    static result_type apply( contiguous_matrix<T,S,column_major> m ) {
      return result_type( m.ptr(), m.num_columns()*m.num_rows() ) ;
    }
  } ;

} // namespace glas2

#endif
