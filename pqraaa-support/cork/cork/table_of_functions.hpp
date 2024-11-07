//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_table_of_functions_hpp
#define cork_table_of_functions_hpp

#include <cork/basis/tabled_functions.hpp>
#include <cork/vector.hpp>
#include <cork/matrix.hpp>

namespace CORK {

  template <typename Z, typename F>
  auto table_of_functions( CORK::vector<Z> const& z_values, CORK::matrix<F> const& functions_values ) {
    return basis::tabled_functions< CORK::vector<Z>, CORK::matrix<F> >( z_values, functions_values ) ;
  }

} // namespace CORK

#endif
