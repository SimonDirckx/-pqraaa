//
// Copyright (c) 2015
// Karl Meerbergen
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
//

#ifndef ARPACK_NO_OP_HPP
#define ARPACK_NO_OP_HPP

#include <cassert>

namespace ARPACK {

  struct no_op {
    template <typename T>
    void operator() ( int n, T* x, T* y ) const {
      assert( false ) ;
    }
  } ;

} // namespace ARPACK

#endif

