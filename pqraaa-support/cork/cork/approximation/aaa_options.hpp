//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_aaa_options_hpp
#define cork_approximation_aaa_options_hpp

#include <cassert>
#include <cmath>
#include <complex>
#include <limits>

namespace CORK { namespace approximation {

  template <typename T=double>
  struct aaa_options {
    typedef T                       value_type ;
    typedef decltype(std::abs(T())) real_type ;

    aaa_options()
    : n_points( 100000 )
    , max_degree( 100 )
    , tolerance( std::numeric_limits<real_type>::epsilon()*100. )
    , function_drop_tol( tolerance )
    , cheap_correction_tolerance (0)
    , check_poles( false )
    , fast( true )
    , lawson_iterations( 0 )
    , debug_level(0)
    {}

    int       n_points ;
    int       max_degree ;
    real_type tolerance ;
    real_type function_drop_tol ;
    real_type cheap_correction_tolerance ;
    bool      check_poles ;
    bool      fast ;
    int       lawson_iterations ;
    int       debug_level ;
  };

} } // namespace CORK::approximation

#endif
