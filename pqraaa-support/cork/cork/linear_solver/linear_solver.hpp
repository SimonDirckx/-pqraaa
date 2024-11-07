//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_linear_solver_linear_solver_hpp
#define cork_linear_solver_linear_solver_hpp

#include <cork/utility/ref.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace linear_solver {

  template <typename BasisValueType, typename CoefficientMatrices, typename EnableIf=void>
  class linear_solver_traits
  {
  } ; // class linear_solver


  template <typename T, typename CoefficientMatrices>
  auto make_linear_solver( CoefficientMatrices const& b ) {
    return linear_solver_traits<T,typename deref_type<CoefficientMatrices>::type>::apply( b ) ;
  }

} } // namespace CORK::coefficient_matrices

#endif
