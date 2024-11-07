//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas_toolbox_tensor_algorithm_fiber_hpp
#define glas_toolbox_tensor_algorithm_fiber_hpp

#include <glas/config.hpp>
#define GLAS_MY_INCLUDE_LEVEL GLAS_INCLUDE_ALGORITHM
#include <glas/concept/check_include_level.hpp>

//#ifdef GLAS_BACKEND_GLAS
//#include <glas/toolbox/tensor/backend/glas/algorithm/fiber.hpp>
//#endif
#include <glas/toolbox/tensor/view/fiber_view.hpp>

#undef GLAS_INCLUDE_LEVEL
#define GLAS_INCLUDE_LEVEL GLAS_INCLUDE_ALGORITHM

namespace glas {

  template <typename T, typename I, typename C>
  fiber_view<T,C> fiber( T& tensor, I index, C const& coordinates ) {
    return fiber_view<T,C>( tensor, index, coordinates ) ;
  }
}

#endif
