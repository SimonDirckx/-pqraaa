//
// Copyright (c) 2015
// Karl Meerbergen
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
//

#ifndef ARPACK_DATA_HPP
#define ARPACK_DATA_HPP

#include <boost/numeric/bindings/value_type.hpp>
#include <boost/numeric/bindings/size.hpp>
#include <boost/numeric/bindings/begin.hpp>
#include <string>
#include <cassert>

namespace ARPACK {

  namespace bindings = boost::numeric::bindings ;

  template <typename T>
  struct data {
    typedef T value_type ;

    template <typename V0, typename V, typename Work>
    data( std::string const& which, int nev, int max_restarts, T const& tol, V0& v0, V& v, Work& work, int info=0 )
    : bmat( '0' )
    , n_( bindings::size(v0) )
    , nev_( nev )
    , tol_( tol )
    , resid_( bindings::begin_value(v0) )
    , ncv_( bindings::size_column(v) )
    , v_( bindings::begin_value(v) )
    , ldv_( bindings::stride_major(v) )
    , workd_( bindings::begin_value(work) )
    , workl_( workd_ + 3 * n_ )
    , lworkl_( bindings::size(work) - 3*n_ )
    , info_( info )
    {
      assert( bindings::size(v0) == bindings::size_row(v) ) ;
      assert( which.length()==2 ) ;
      which_[0] = which[0] ; which_[1] = which[1] ;
      iparam_[3] = 1 ;
      iparam_[2] = max_restarts ;
    }

    // Data
    char bmat ;
    int  n_ ;
    char which_[2] ;
    int  nev_ ;
    T    tol_ ;
    T*   resid_ ;
    int  ncv_ ;
    T*   v_ ;
    int  ldv_ ;
    int  iparam_[11] ;
    int  ipntr_[14] ;
    T*   workd_ ;
    T*   workl_ ;
    int  lworkl_ ;
    int  info_ ;
  } ;

} // namespace ARPACK

#endif

