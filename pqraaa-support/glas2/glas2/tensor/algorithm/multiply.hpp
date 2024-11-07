//  (C) Copyright Karl Meerbergen 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_tensor_algorithm_multiply_hpp
#define glas2_tensor_algorithm_multiply_hpp

#include <glas2/tensor/expression/tensor_matrix_multiply.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/tensor/concept/dense_tensor.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>

namespace glas2 {

  template <typename T, typename M, typename S>
  typename std::enable_if< is<DenseTensor,T>::value && is<DenseMatrix,M>::value && std::is_integral<S>::value
                         , tensor_matrix_multiply<T,M>
                         >::type multiply( T const& t, M const& m, S mode ) {
    return tensor_matrix_multiply<T,M>( t, m, mode ) ;
  }

} // namespace glas2

#endif
