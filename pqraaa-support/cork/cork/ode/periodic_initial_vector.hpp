//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_ode_constant_initial_vector_hpp
#define cork_ode_constant_initial_vector_hpp

#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>

namespace CORK { namespace ode {

  template <typename X0, typename Basis, typename Shift, typename V0>
  void periodic_initial_vector( X0 const& x0, Basis const& basis4cork, Shift const& shift, V0& v0 ) {
    auto handle = basis4cork.template handle<typename V0::value_type>() ;
    typedef typename decltype(handle)::value_type value_type ;
    handle.shift( shift ) ;
    glas2::vector< value_type > rhs( 1 ) ; rhs(0) = -1.0 ;
    glas2::matrix< value_type > Z( basis4cork.size()-1, 1 ) ;
    handle.lower_solve_right_hand_side( rhs, Z ) ;
    handle.solve( Z ) ;
    v0 ( glas2::range(0, x0.size()) ) = x0 ;
    reshape( v0( glas2::range_from_end(x0.size(),0) ), x0.size(), Z.num_rows(), glas2::column_major() ) = outer_prod( x0, Z( glas2::all(), 0 ) ) ;
  } // periodic_initial_vector()
   
} } // namespace CORK::ode


#endif
