//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_algorithm_transpose_hpp
#define glas2_matrix_algorithm_transpose_hpp

#include <glas2/matrix/type/contiguous_matrix.hpp>
#include <glas2/matrix/type/strided_matrix.hpp>
#include <glas2/matrix/type/double_strided_matrix.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/matrix/concept/matrix_transform.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>

namespace glas2 {

  struct transpose_tag {} ;

  template <typename M>
  typename std::enable_if< is<DenseMatrix,M>::value, typename matrix_transform< transpose_tag, M >::result_type >::type transpose( M m ) {
    return matrix_transform< transpose_tag, M>::apply(m) ;
  }

  template <typename T, typename S>
  struct matrix_transform< transpose_tag, contiguous_matrix<T,S,column_major> > {
    typedef contiguous_matrix<T,S,row_major> result_type ;

    static result_type apply( contiguous_matrix<T,S,column_major> m ) {
      return result_type( m.ptr(), m.num_columns(), m.num_rows() ) ;
    }
  } ;

  template <typename T, typename S>
  struct matrix_transform< transpose_tag, contiguous_matrix<T,S,row_major> > {
    typedef contiguous_matrix<T,S,column_major> result_type ;

    static result_type apply( contiguous_matrix<T,S,row_major> m ) {
      return result_type( m.ptr(), m.num_columns(), m.num_rows() ) ;
    }
  } ;

  template <typename T, typename S>
  struct matrix_transform< transpose_tag, strided_matrix<T,S,column_major> > {
    typedef strided_matrix<T,S,row_major> result_type ;

    static result_type apply( strided_matrix<T,S,column_major> m ) {
      return result_type( m.ptr(), m.stride(), m.num_columns(), m.num_rows() ) ;
    }
  } ;

  template <typename T, typename S>
  struct matrix_transform< transpose_tag, strided_matrix<T,S,row_major> > {
    typedef strided_matrix<T,S,column_major> result_type ;

    static result_type apply( strided_matrix<T,S,row_major> m ) {
      return result_type( m.ptr(), m.stride(), m.num_columns(), m.num_rows() ) ;
    }
  } ;

  template <typename T, typename S>
  struct matrix_transform< transpose_tag, double_strided_matrix<T,S,column_major> > {
    typedef double_strided_matrix<T,S,row_major> result_type ;

    static result_type apply( double_strided_matrix<T,S,column_major> m ) {
      return result_type( m.ptr(), m.stride_columns(), m.stride_rows(), m.num_columns(), m.num_rows() ) ;
    }
  } ;

  template <typename T, typename S>
  struct matrix_transform< transpose_tag, double_strided_matrix<T,S,row_major> > {
    typedef double_strided_matrix<T,S,column_major> result_type ;

    static result_type apply( double_strided_matrix<T,S,row_major> m ) {
      return result_type( m.ptr(), m.stride_columns(), m.stride_rows(), m.num_columns(), m.num_rows() ) ;
    }
  } ;

} // namespace glas2

#endif
