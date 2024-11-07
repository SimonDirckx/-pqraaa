//  (C) Copyright Karl Meerbergen & Dries De Samblanx, 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_scalar_functions_hpp
#define cork_scalar_functions_hpp

#include <cork/basis/set_of_functions.hpp>
#include <cork/basis/tabled_functions.hpp>
#include <cassert>

namespace CORK {

  template <typename Domain, typename Functions>
  auto scalar_functions( Domain const& domain, Functions const& F, int m ) {
    return basis::set_of_functions<Domain,Functions>( domain, F, m ) ;
  } // scalar_functions

  template <typename Domain, typename Table>
  auto scalar_functions( Domain const& domain, Table const& F ) {
    return basis::set_of_functions<Domain,Functions>( domain, F, m ) ;
  } // scalar_functions

} // namespace CORK

#endif

#endif
