//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_tensor_algorithm_cast_to_dense_matrix_hpp
#define glas2_tensor_algorithm_cast_to_dense_matrix_hpp

#include <glas2/tensor/concept/contiguous_dense_tensor.hpp>
#include <glas2/matrix/type/contiguous_matrix.hpp>
#include <glas2/matrix/type/double_strided_matrix.hpp>
#include <glas2/matrix/type/strided_matrix.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename T>
  typename std::enable_if< is<DenseTensor,T>::value, strided_matrix<> >::type
                         >::type cast_to_dense_matrix( dense_unfolding<T,S> unfolding ) {
    return double_strided_matrix<typename T::value_type
                                , typename T::size_type
                                , column_major
                                >( tensor.ptr(), tensor.stride()(0), tensor.stride()(1), tensor.shape()(0), tensor.shape()(1) ) ;
  }

} // namespace glas2

#endif
