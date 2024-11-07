//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_sparse_nonlinear_hpp
#define cork_matrix_valued_function_sparse_nonlinear_hpp

#include <cork/matrix_valued_function/nonlinear.hpp>
#include <cork/sparse.hpp>
#include <cassert>
#include <type_traits>
#include <vector>

namespace CORK { namespace matrix_valued_function {

//  template <typename Basis, typename FunctionSequence, typename CoefficientMatricesSequence>
//  decltype (auto) make_sparse_nonlinear_lvalue( Basis const& basis, FunctionSequence const& function_sequence, CoefficientMatricesSequence const& sequence )  {
//    coefficient_matrices::sparse< CoefficientMatricesSequence > matrices( sequence ) ;
//    return make_nonlinear_lvalue( basis, function_sequence, matrices ) ;
//  }

  template <typename Basis, typename FunctionSequence, typename CoefficientMatricesSequence>
  decltype (auto) make_sparse_nonlinear( Basis const& basis, FunctionSequence const& function_sequence, CoefficientMatricesSequence const& sequence )  {
    coefficient_matrices::sparse< CoefficientMatricesSequence > matrices( sequence ) ;
    basis::union_of_functions<Basis, FunctionSequence> polybasis(basis, function_sequence);
    return make_matrix_polynomial( polybasis, matrices )
//    return make_nonlinear( basis, function_sequence, matrices ) ;
  }

} } // namespace CORK::matrix_valued_function

#endif
