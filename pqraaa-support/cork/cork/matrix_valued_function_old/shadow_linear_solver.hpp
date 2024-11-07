//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_shadow_linear_solver_hpp
#define cork_matrix_valued_function_shadow_linear_solver_hpp

#include <cork/matrix_valued_function/matrix_polynomial_with_shadow.hpp>
#include <cork/matrix_valued_function/linear_solver.hpp>

namespace CORK { namespace matrix_valued_function {

  template <typename ValueType, typename MatrixPolynomial, typename Shadow, typename CoefficientMatrices>
  struct linear_solver_struct< ValueType, matrix_polynomial_with_shadow<MatrixPolynomial,Shadow>, CoefficientMatrices > {
    typedef linear_solver_struct< ValueType, Shadow, typename Shadow::coefficient_matrices_type > struct_type ;
    typedef typename struct_type::value_type   value_type ;
    typedef typename struct_type::solver_type  solver_type ;
    typedef typename struct_type::type         type ;

    static type apply( matrix_polynomial_with_shadow<MatrixPolynomial,Shadow> const& mp ) {
      return struct_type::apply( mp.shadow() ) ;
    }
  } ; // struct linear_solver_struct

} } // namespace CORK::matrix_valued_function

#endif
