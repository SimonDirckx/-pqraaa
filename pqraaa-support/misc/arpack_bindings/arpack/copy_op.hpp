//
// Copyright (c) 2015
// Karl Meerbergen
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
//

#ifndef ARPACK_COPY_OP_HPP
#define ARPACK_COPY_OP_HPP

#include <cassert>
#include <algorithm>

namespace ARPACK {

  struct copy_op {
    template <typename T>
    void operator() ( int n, T* x, T* y ) const {
      std::copy( x, x+n, y ) ;
    }

    static int const info = 0 ;
  } ;

} // namespace ARPACK

#endif

