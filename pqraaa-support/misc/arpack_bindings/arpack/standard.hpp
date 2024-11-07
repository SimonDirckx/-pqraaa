//
// Copyright (c) 2015
// Karl Meerbergen
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
//

#ifndef ARPACK_STANDARD_HPP
#define ARPACK_STANDARD_HPP

#include <arpack/no_op.hpp>
#include <string>
#include <cassert>

namespace ARPACK {

  //
  // Standard eigenvalue problem
  //   A * x = lambda * x
  // solved by Krylov space with A
  //
  // op(x,y) computes y = A * x
  //
  template <typename T, typename Op>
  struct standard {
    standard( Op const& op )
    : op_( op )
    {}

    static const char bmat = 'I' ;
    static const int  mode = 1 ;

    T sigma() const {return T() ;}
    T sigma_r() const {return T() ;}
    T sigma_i() const {return T() ;}

    Op const& op_ ;
    no_op matvec_b_ ;
    no_op initial_op_ ;
  } ;

  template <typename T, typename Op>
  standard<T,Op> make_standard( T, Op const& op ) {
    return standard<T,Op>( op ) ;
  }


} // namespace ARPACK

#endif

