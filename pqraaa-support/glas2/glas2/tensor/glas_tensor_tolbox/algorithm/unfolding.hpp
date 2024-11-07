//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas_toolbox_tensor_algorithm_unfolding_hpp
#define glas_toolbox_tensor_algorithm_unfolding_hpp

#include <glas/config.hpp>
#define GLAS_MY_INCLUDE_LEVEL GLAS_INCLUDE_ALGORITHM
#include <glas/concept/check_include_level.hpp>

#undef GLAS_INCLUDE_LEVEL
#define GLAS_INCLUDE_LEVEL GLAS_INCLUDE_ALGORITHM

namespace glas {

  template <typename T, typename R, typename C>
  unfolding_view<T,R,C> unfolding( T& tensor, R const& rows, C const& columns ) {
    return unfolding_view<T,R,C>( tensor, rows, columns ) ;
  }
}

#endif
