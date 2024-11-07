//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_set_of_functions_hpp
#define cork_set_of_functions_hpp

#include <cork/basis/set_of_functions.hpp>
#include <cork/vector.hpp>
#include <cork/matrix.hpp>

namespace CORK {

  template <typename ValueType, typename F>
  auto set_of_functions( F const& fun ) {
    return basis::set_of_functions< ValueType, F >( fun, 0 ) ;
  }

} // namespace CORK

#endif
