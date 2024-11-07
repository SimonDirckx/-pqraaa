//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_tensor_algorithm_tensor2matrix_hpp
#define glas2_tensor_algorithm_tensor2matrix_hpp

#include <glas2/tensor/concept/contiguous_dense_tensor.hpp>
#include <glas2/matrix/type/contiguous_matrix.hpp>
#include <glas2/matrix/type/double_strided_matrix.hpp>
#include <glas2/matrix/type/strided_matrix.hpp>
#include <glas2/tensor/type/dense_unfolding.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename T>
  typename std::enable_if< is<ContiguousDenseTensor,T>::value
                         , contiguous_matrix<typename T::value_type, typename T::size_type, column_major >
                         >::type tensor2matrix( T tensor ) {
    assert( 2==tensor.order() ) ;
    return contiguous_matrix<typename T::value_type, typename T::size_type, column_major >( tensor.ptr(), tensor.shape()(0), tensor.shape()(1) ) ;
  }

  template <typename T>
  typename std::enable_if< is<StridedDenseTensor,T>::value
                         , double_strided_matrix<typename T::value_type, typename T::size_type, column_major >
                         >::type tensor2matrix( T tensor ) {
    assert( 2==tensor.order() ) ;
    return double_strided_matrix<typename T::value_type
                                , typename T::size_type
                                , column_major
                                >( tensor.ptr(), tensor.stride()(0), tensor.stride()(1), tensor.shape()(0), tensor.shape()(1) ) ;
  }

  template <class T, class S>
  double_strided_matrix<T, S, column_major > tensor2matrix( dense_unfolding<T,S> unfolding ) {
#ifndef NDEBUG
    {
      const int stride = unfolding.begin_column(1)-unfolding.begin_column(0) ;
      for (int i=1; i<unfolding.num_columns()-1; ++i) {
        assert( stride==unfolding.begin_column(i+1)-unfolding.begin_column(i) ) ;
      }
    }
#endif
    return double_strided_matrix< T, S, column_major >( unfolding.ptr(), unfolding.stride_rows(), unfolding.begin_column(1)-unfolding.begin_column(0)
                                 , unfolding.num_rows(), unfolding.num_columns()
                                 ) ;
  }

} // namespace glas2

#endif
