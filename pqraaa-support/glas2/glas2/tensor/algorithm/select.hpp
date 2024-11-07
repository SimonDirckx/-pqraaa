//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_tensor_algorithm_select_hpp
#define glas2_tensor_algorithm_select_hpp

#include <glas2/vector/type/all.hpp>
#include <glas2/vector/type/slice.hpp>
#include <glas2/tensor/concept/dense_tensor.hpp>
#include <glas2/tensor/view/tensor_selection.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>

namespace glas2 {

  template <typename T, typename S, typename Mode>
  typename std::enable_if< is<DenseTensor,T>::value, typename tensor_selection< T, S, Mode >::result_type >::type select( T t, S s, Mode m ) {
    return tensor_selection< T,S,Mode >::apply(t,s,m) ;
  }

} // namespace glas2

#endif
