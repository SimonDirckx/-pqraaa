//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_coefficient_matrices_linear_solver_hpp
#define cork_coefficient_matrices_linear_solver_hpp

#include <cassert>
#include <type_traits>

namespace CORK { namespace coefficient_matrices {

  template <typename BasisValueType, typename CoefficientMatrices>
  class linear_solver
  {
  } ; // class linear solver


  template <typename T, typename CoefficientMatrices>
  decltype(auto) make_linear_solver( CoefficientMatrices const& b ) {
    return linear_solver<T,typename std::decay<CoefficientMatrices>::type >( b ) ;
  }

} } // namespace CORK::coefficient_matrices

#endif
