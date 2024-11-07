//  (C) Copyright Karl Meerbergen 2022.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_matrix_valued_function_with_solver_hpp
#define cork_matrix_valued_function_matrix_valued_function_with_solver_hpp

#include <cork/options/value_of.hpp>
#include <cork/options/relative_shift_modifier.hpp>
#include <cork/options/absolute_shift_modifier.hpp>
#include <cork/exception/linear_solver_failure.hpp>
#include <cork/linear_solver/linear_solver.hpp>
#include <cork/matrix_valued_function/matrix_polynomial.hpp>
#include <cork/matrix_valued_function/nonlinear.hpp>
#include <cork/basis4cork/basis4cork.hpp>
#include <cork/utility/ref.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <type_traits>
#include <cassert>
#include <typeinfo>

namespace CORK { namespace matrix_valued_function {

  template <typename MatVec, typename Solver>
  class matrix_valued_function_with_solver
  {
    public:
      typedef typename CORK::deref_type< CoefficientMatricesLinearSolver >::type linear_solver_type ;
      typedef typename CORK::deref_type< CoefficientMatrices >::type             coefficient_matrices_type ;
      typedef typename CORK::deref_type< Basis >::type                           sequence_of_functions_type ;
      typedef typename linear_solver_type::size_type                             size_type ;

    public:
      typedef typename linear_solver_type::value_type value_type ;

    public:
      explicit matrix_valued_function_with_solver( MatVec const& matvec, Solver solver )
      : matvec_( matvec )
      , solver_( solver )
      {}

    private:
      size_type num_terms() const { return num_terms_ ; }

    public:
      size_type size() const { return size_ ; }

    private:
      MatVec matvec_ ;
      Solver solver_ ;
  } ; // matrix_valued_function_with_solver

} } // namespace CORK::matrix_valued_function

#endif
