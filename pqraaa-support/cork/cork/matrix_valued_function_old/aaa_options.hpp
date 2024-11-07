//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_aaa_options_hpp
#define cork_matrix_valued_function_aaa_options_hpp

#include <cork/approximation/aaa_options.hpp>

namespace CORK { namespace matrix_valued_function {

  template <typename T=double>
  struct aaa_options
  : approximation::aaa_options<T>
  {
    aaa_options()
    : test_error( false )
    , compute_poles( false )
    {}

    bool test_error ;
    bool compute_poles ;
  };

} } // namespace CORK::matrix_valued_function

#endif
