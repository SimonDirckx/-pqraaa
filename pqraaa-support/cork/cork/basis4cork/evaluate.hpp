//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_basis4cork_evaluate_hpp
#define cork_basis4cork_evaluate_hpp

#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <string>
#ifndef NDEBUG
#include <limits>
#endif

namespace CORK { namespace basis4CORK {

  // Only valid for bases for which the size is equal to num_terms
  template <typename Basis, typename Coefs>
  void evaluate( Basis const& basis, Coefs& coefs ) {
    glas2::matrix< typename Coefs::value_type > z( coefs.size(), 1 ) ;
    z(0,0) = 1.0 ;
    auto zr = z(glas2::range_from_end(1,0), glas2::all()) ;
    basis.lower_solve_right_hand_side(z(0,glas2::all()),zr) ;
    basis.solve(zr) ;
    zr *= -1.0 ;
    coefs = z(glas2::all(),0) ;
  } // evaluate()

} } // namespace CORK::basis4cork

#endif
