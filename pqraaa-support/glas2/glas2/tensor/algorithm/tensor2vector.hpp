//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_tensor_algorithm_tensor2vector_hpp
#define glas2_tensor_algorithm_tensor2vector_hpp

#include <glas2/tensor/concept/contiguous_dense_tensor.hpp>
#include <glas2/tensor/concept/strided_dense_tensor.hpp>
#include <glas2/vector/type/contiguous_vector.hpp>
#include <glas2/vector/type/double_strided_vector.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename T>
  typename std::enable_if< is<ContiguousDenseTensor,T>::value
                         , contiguous_vector<typename T::value_type, typename T::size_type >
                         >::type tensor2vector( T tensor ) {
    assert( 1==tensor.order() ) ;
    return contiguous_vector<typename T::value_type, typename T::size_type >( tensor.ptr(), tensor.shape()(0) ) ;
  }

  template <typename T>
  typename std::enable_if< is<StridedDenseTensor,T>::value
                         , double_strided_vector<typename T::value_type, typename T::size_type, column_major >
                         >::type tensor2vector( T tensor ) {
    assert( 1==tensor.order() ) ;
    return strided_vector<typename T::value_type
                         , typename T::size_type
                         >( tensor.ptr(), tensor.shape()(0), tensor.stride()(0) ) ;
  }

} // namespace glas2

#endif
