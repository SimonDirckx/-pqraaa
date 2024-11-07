//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_linearization_clone_cork_linearization_hpp
#define cork_linearization_clone_cork_linearization_hpp

#include <cork/linearization/cork_linearization.hpp>
#include <cork/matrix_iterator/matrix_iterator.hpp>
#include <cork/basis4cork/basis4cork.hpp>

namespace CORK { namespace linearization {

  template <typename CORKLinearization>
  decltype(auto) clone_cork_linearization( CORKLinearization& linearization ) {
    auto linear_solver = linearization.linear_solver().clone() ;
    return cork_linearization< decltype(linearization.basis()) const&
                             , decltype(linearization.matrix_iterator())&
                             , decltype(linearization.coefficient_matrices()) const&
                             , decltype(linear_solver)
                             >( linearization.basis(), linearization.matrix_iterator(), linearization.coefficient_matrices(), linear_solver, linearization.information() ) ;
  }

} } // namespace CORK::linearization

#endif
