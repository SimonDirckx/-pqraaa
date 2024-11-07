//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_potential_theory_options_hpp
#define cork_matrix_valued_function_potential_theory_options_hpp

#include <cork/approximation/potential_theory_options.hpp>

namespace CORK { namespace matrix_valued_function {

  template <typename T=double>
  struct potential_theory_options
  : approximation::potential_theory_options<T>
  {
    potential_theory_options()
    : test_error( false )
    {}

    bool test_error ;
  };

} } // namespace CORK::matrix_valued_function

#endif
