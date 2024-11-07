//
// Copyright (c) 2015
// Karl Meerbergen
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
//

#ifndef ARPACK_RANDOM_VECTOR_OP_HPP
#define ARPACK_RANDOM_VECTOR_OP_HPP

#include <cassert>
#include <algorithm>

namespace ARPACK {

  struct random_vector_op {
    template <typename T>
    void operator() ( int n, T* x, T* y ) const {
      assert( false ) ;
    }

    static int const info = 1 ;
  } ;

} // namespace ARPACK

#endif

