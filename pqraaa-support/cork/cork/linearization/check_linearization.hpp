//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_linearization_check_linearization_hpp
#define cork_linearization_check_linearization_hpp

#include <cmath>
#include <type_traits>

namespace CORK { namespace linearization {

  //
  // Computes ||P(shift)* x - lin * (Phi\otimes x)||
  //
  template <typename Linearization, typename Shift, typename X>
  decltype(std::abs(Shift())) check_linearization( Linearization const& lin, Shift const& shift, X const& x ) {
    typedef typename std::common_type<Shift,typename X::value_type>::type value_type ;

    // Vector Phi
    auto basis_handle = lin.basis().template handle<value_type>() ;
    basis_handle.shift( shift ) ;
    glas2::matrix< value_type > coefs( lin.size_of_basis(), 1 ) ;
    coefs(0, 0) = 1. ;
    basis_handle.lower_solve_right_hand_side( -coefs(0, glas2::all()), coefs( glas2::range_from_end(1,0), glas2::all()) ) ;
    basis_handle.solve( coefs( glas2::range_from_end(1,0), glas2::all()) ) ;
    //std::cout << "coefs " << coefs << std::endl ;

    // Vector Phi * x
    glas2::vector< value_type > xx( x.size() * lin.size_of_basis() ) ;
    xx = glas2::kron( coefs(glas2::all(),0), x ) ;

    auto multiply_handle = lin.template multiply_handle<value_type>() ;
    glas2::vector< value_type > yy( xx.size() ) ;
    glas2::vector< value_type > yy_temp( xx.size() ) ;
    multiply_handle.multiply_A( reshape(xx,lin.size_of_basis(),x.size(),glas2::row_major()), reshape(yy,lin.size_of_basis(),x.size(),glas2::row_major()) ) ;
    multiply_handle.multiply_B( reshape(xx,lin.size_of_basis(),x.size(),glas2::row_major()), reshape(yy_temp,lin.size_of_basis(),x.size(),glas2::row_major()) ) ;
    yy -= shift * yy_temp ;
    //std::cout << "yy " << yy << std::endl ;

    yy(glas2::range(0,x.size())) *= -1.0 ;
    lin.matrix_polynomial().multiply_add( shift, x, yy(glas2::range(0,x.size())) ) ;
    //std::cout << "yy " << reshape(yy,lin.size_of_basis(),x.size(),glas2::row_major()) << std::endl ;

    return norm_2( yy ) ;
  } // check_linearization
   
} } // namespace CORK::linearization


#endif
