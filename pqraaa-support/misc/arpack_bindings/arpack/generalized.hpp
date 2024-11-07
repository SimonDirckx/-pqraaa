//
// Copyright (c) 2015
// Karl Meerbergen
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
//

#ifndef ARPACK_GENERALIZED_HPP
#define ARPACK_GENERALIZED_HPP

#include <string>
#include <cassert>

namespace ARPACK {

  //
  // Generalized eigenvalue problem
  //   A * x = lambda * B * x
  // B is assumed positive definite.
  //
  // The problem is solved by Krylov space with B^{-1} * A
  //
  // op(x,y) computes y = B^{-1} * A * x
  // matvec_B(x,y) computes y = B * x
  //
  template <typename T, typename Op, typename OpB>
  struct generalized {
    generalized( Op const& op, OpB const& matvec_B )
    : op_( op )
    , matvec_B_( matvec_B )
    {}

    static const char bmat = 'G' ;
    static const int  mode = 2 ;

    T  sigma() { return T() ;}
    T sigma_r() const {return T() ;}
    T sigma_i() const {return T() ;}

    Op const&  op_ ;
    OpB cosnt& matvec_B_ ;
  } ;


  template <typename T, typename Op, typename OpB>
  generalized<T,Op,OpB> make_generalized( T, Op const& op, OpB const& matvec_B ) {
    return generalized<T,Op,OpB>( op, matvec_B ) ;
  }

} // namespace ARPACK

#endif

