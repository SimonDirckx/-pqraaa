//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_algorithm_diagonal_hpp
#define glas2_matrix_algorithm_diagonal_hpp

#include <glas2/vector/type/strided_vector.hpp>
#include <glas2/matrix/type/contiguous_matrix.hpp>
#include <glas2/matrix/type/strided_matrix.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/matrix/concept/matrix_transform.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cmath>

namespace glas2 {

  struct diagonal_tag {} ;

  template <typename M>
  typename std::enable_if< is<DenseMatrix,M>::value, typename matrix_transform< diagonal_tag, M >::result_type >::type diagonal( M m, typename M::size_type d=0 ) {
    return matrix_transform< diagonal_tag, M>::apply(m,d) ;
  }

  template <typename T, typename S>
  struct matrix_transform< diagonal_tag, contiguous_matrix<T,S,column_major> > {
    typedef strided_vector<T,S> result_type ;

    static result_type apply( contiguous_matrix<T,S,column_major> m, S d ) {
      if (d>=0)
        return result_type( &m(0,d), std::min(m.num_columns()-d,m.num_rows()), m.num_rows()+1 ) ;
      else
        return result_type( &m(-d,0), std::min(m.num_columns(),m.num_rows()+d), m.num_rows()+1 ) ;
    }
  } ;

  template <typename T, typename S>
  struct matrix_transform< diagonal_tag, contiguous_matrix<T,S,row_major> > {
    typedef strided_vector<T,S> result_type ;

    static result_type apply( contiguous_matrix<T,S,row_major> m, S d ) {
      if (d>=0)
        return result_type( &m(0,d), std::min(m.num_columns()-d,m.num_rows()), m.num_columns()+1 ) ;
      else
        return result_type( &m(-d,0), std::min(m.num_columns(),m.num_rows()+d), m.num_columns()+1 ) ;
    }
  } ;

  template <typename T, typename S>
  struct matrix_transform< diagonal_tag, strided_matrix<T,S,column_major> > {
    typedef strided_vector<T,S> result_type ;

    static result_type apply( strided_matrix<T,S,column_major> m, S d ) {
      if (d>=0)
        return result_type( &m(0,d), std::min(m.num_columns()-d,m.num_rows()), m.stride()+1 ) ;
      else
        return result_type( &m(-d,0), std::min(m.num_columns(),m.num_rows()+d), m.stride()+1 ) ;
    }
  } ;

  template <typename T, typename S>
  struct matrix_transform< diagonal_tag, strided_matrix<T,S,row_major> > {
    typedef strided_vector<T,S> result_type ;

    static result_type apply( strided_matrix<T,S,row_major> m, S d ) {
      if (d>=0)
        return result_type( &m(0,d), std::min(m.num_columns()-d,m.num_rows()), m.stride()+1 ) ;
      else
        return result_type( &m(-d,0), std::min(m.num_columns(),m.num_rows()+d), m.stride()+1 ) ;
    }
  } ;

} // namespace glas2

#endif
