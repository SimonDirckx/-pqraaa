//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_basis4cork_solve_2_by_2_hpp
#define cork_basis4cork_solve_2_by_2_hpp

#include <cork/basis4cork/basis4cork.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <glas2/bindings/matrix.hpp>
#include <glas2/bindings/vector.hpp>
#include <cassert>
#include <string>
#ifndef NDEBUG
#include <limits>
#endif

namespace CORK { namespace basis4CORK {

  template <typename LU, typename Z>
  void solve_2_by_2( LU& lu, Z z ) {
    assert( lu.num_rows()==2 ) ;
    assert( lu.num_columns()==2 ) ;
    bool pivot = false ;
    if (std::abs(lu(0,0))<std::abs(lu(1,0))) pivot = true ;

    if (pivot) {
      glas2::swap( lu(0,glas2::all()), lu(1,glas2::all()) ) ;
      glas2::swap( z(0,glas2::all()), z(1,glas2::all()) ) ;
    }

    // Factorize
    lu(1,0) /= lu(0,0) ;
    lu(1,1) -= lu(1,0) * lu(0,1) ;

    // Forward solve
    z(1,glas2::all()) -= lu(1,0) * z(0,glas2::all()) ;

    // Backward solve
    z(1,glas2::all()) /= lu(1,1) ;
    z(0,glas2::all()) -= lu(0,1) * z(1,glas2::all()) ;
    z(0,glas2::all()) /= lu(0,0) ;
  } ; // solve_2_by_2()

} } // namespace CORK::basis4cork

#endif
