//  (C) Copyright Karl Meerbergen 2011. 
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_iterative_krylov_sequential_hpp
#define glas2_iterative_krylov_sequential_hpp

#include <glas2/vector.hpp>
#include <cmath>

namespace glas2 { namespace iterative {


  struct sequential {
    template <typename X, typename Y>
    decltype( typename X::value_type() + typename Y::value_type() ) inner_prod( X x, Y y ) const {
        return glas2::inner_prod( x, y ) ;
    }

    template <typename X>
    decltype( std::abs( typename X::value_type() ) ) norm_2( X x ) const {
        return glas2::norm_2( x ) ;
    }
  } ;

} }

#endif
