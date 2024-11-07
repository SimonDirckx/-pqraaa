//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_USE_SHADOWED_NONLINEAR_MATRIX_HPP
#define CORK_OPTIONS_USE_SHADOWED_NONLINEAR_MATRIX_HPP

#include<limits>
#include<tuple>

namespace CORK { namespace options {

  class use_shadowed_nonlinear_matrix {
    public:
      use_shadowed_nonlinear_matrix()
      {}

      template <typename ...Ts>
      use_shadowed_nonlinear_matrix( std::tuple<Ts...> const& options )
      {}
  } ;

} } // CORK::options

#endif
