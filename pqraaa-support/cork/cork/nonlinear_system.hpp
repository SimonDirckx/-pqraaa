//  (C) Copyright Karl Meerbergen & Dries De Samblanx, 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_nonlinear_system_hpp
#define cork_nonlinear_system_hpp

#include <cork/basis/monomial.hpp>
#include <cork/set_of_functions.hpp>
#include <cork/basis/tabled_functions.hpp>
#include <cork/basis/union_of_functions.hpp>
#include <cork/concept/multiply_add.hpp>
#include <cork/concept/set_of_scalar_functions.hpp>
#include <cork/concept/scalar_functions_sequence.hpp>
#include <cork/concept/solve.hpp>
#include <cork/matrix_valued_function/matrix_polynomial.hpp>
#include <cork/coefficient_matrices/matrices_by_functions.hpp>
#include <cork/user_defined/linear_solver.hpp>
#include <cork/basis/functions.hpp>
#include <cork/basis/union_of_bases.hpp>
#include <cassert>
#include <typeinfo>

namespace CORK {

 /**
  * Returns a nonlinear_system object that represents the nonlinear eigenvalue problem.
  *
  * @param n                  The size of the problem.
  * @param degree             The degree of the polynomial part.
  * @param scalar_function    The scalar functions of the non-polynomial part.
  *                           Can be given as a sequence of functions, functors or lambdas or as a table of function values.
  * @param MultiplyAdd        A function that satisfies the interface multiply_add(i, alpha, n, x, y);
  * @param solve              A function that satisfies the interface solve(shift, n, x, is_new_shift);
  */

  template <typename ValueType, typename ScalarFunctions, typename FunctionsValueType, typename MultiplyAdd, typename Solve>
#ifdef CORK_USE_CONCEPTS
        requires CORK::Arithmetic<ValueType> && CORK::Solve<Solve, ValueType>
                && CORK::MultiplyAdd<MultiplyAdd,ValueType>
                && CORK::SetOfScalarFunctions< ScalarFunctions, FunctionsValueType >
#endif
  auto nonlinear_system( int n, int degree, int n_fun, basis::set_of_functions< FunctionsValueType, ScalarFunctions > const& scalar_functions, MultiplyAdd const& multiply_add, Solve& solve ) {
    assert(n>0) ;
    assert(degree>0) ;
    assert(n_fun>0) ;

    basis::set_of_functions< FunctionsValueType, ScalarFunctions > functions( scalar_functions.function(), n_fun ) ;
    CORK::basis::monomial<>                                                            basis       ( degree ) ;
    CORK::basis::union_of_functions<decltype(basis), decltype(functions)>              polybasis   ( basis, functions );
    coefficient_matrices::matrices_by_functions< ValueType, MultiplyAdd, Solve >       coefs       ( n, degree+functions.num_terms(), multiply_add, solve ) ;

    return CORK::matrix_valued_function::make_matrix_polynomial( polybasis, coefs );
  } // nonlinear_system()

  template <typename ValueType, typename ScalarFunctions, typename MultiplyAdd, typename Solve>
#ifdef CORK_USE_CONCEPTS
        requires CORK::Arithmetic<ValueType> && CORK::Solve<Solve, ValueType>
                && CORK::MultiplyAdd<MultiplyAdd, ValueType>
                && CORK::SetOfScalarFunctions< ScalarFunctions, ValueType >
        //&& CORK::ScalarFunctionsSequence< ScalarFunctions, ValueType >
#endif
  auto nonlinear_system( int n, int degree, int n_fun, ScalarFunctions const& scalar_functions, MultiplyAdd const& multiply_add, Solve& solve
                       , typename std::enable_if< !CORK::basis::is_tabled_functions<ScalarFunctions>::value, int>::type i=0 ) {
    assert(n>0) ;
    assert(degree>0) ;
    assert(n_fun>0) ;

    CORK::basis::monomial<>                                                            basis       ( degree ) ;
    CORK::basis::set_of_functions< ValueType, ScalarFunctions >                        functions   ( scalar_functions, n_fun ) ;
    CORK::basis::union_of_functions<decltype(basis), decltype(functions)>              polybasis   ( basis, functions );
    coefficient_matrices::matrices_by_functions< ValueType, MultiplyAdd, Solve >       coefs       ( n, degree+functions.num_terms(), multiply_add, solve ) ;

    return CORK::matrix_valued_function::make_matrix_polynomial( polybasis, coefs );
  } // nonlinear_system()

  template <typename ValueType, typename ScalarFunctions, typename MultiplyAdd, typename Solve>
#ifdef CORK_USE_CONCEPTS
        requires CORK::Arithmetic<ValueType> && CORK::Solve<Solve, ValueType>
                && CORK::MultiplyAdd<MultiplyAdd,ValueType>
#endif
  auto nonlinear_system( int n, int degree, ScalarFunctions const& scalar_functions, MultiplyAdd const& multiply_add, Solve& solve
                                    , typename std::enable_if< CORK::basis::is_tabled_functions<ScalarFunctions>::value, int>::type i=0 ) {
    assert(n>0) ;
    assert(degree>0) ;

    CORK::basis::monomial<>                                                            basis       ( degree ) ;
    CORK::basis::union_of_functions<decltype(basis), decltype(scalar_functions)>       polybasis   ( basis, scalar_functions );
    coefficient_matrices::matrices_by_functions< ValueType, MultiplyAdd, Solve >       coefs       ( n, degree+scalar_functions.num_terms(), multiply_add, solve ) ;

    return CORK::matrix_valued_function::make_matrix_polynomial( polybasis, coefs );
  } // nonlinear_system()

} // namespace CORK

// Insert for dense or sparse matrices if included
#if defined(cork_dense_hpp) || defined(cork_sparse_hpp)

namespace CORK {

  template <typename ValueType, typename MatrixSequence, typename ScalarFunctions>
#ifdef CORK_USE_CONCEPTS
        requires CORK::SetOfScalarFunctions< ScalarFunctions, ValueType >
#endif
  auto nonlinear_system( MatrixSequence const& matrix_sequence, int n_fun, ScalarFunctions const& scalar_functions
                       , typename std::enable_if< !CORK::basis::is_tabled_functions<ScalarFunctions>::value, int>::type i=0 ) {
    assert(n_fun>0) ;
    typedef ValueType value_type ;

    CORK::basis::monomial<>                                                basis       ( matrix_sequence.size()-n_fun-1 ) ;
    CORK::basis::set_of_functions< value_type, ScalarFunctions >           functions   ( scalar_functions, n_fun ) ;
    CORK::basis::union_of_functions<decltype(basis), decltype(functions)>  polybasis   ( basis, functions );
    coefficient_matrices::any_glas< MatrixSequence const& >                coefs       ( matrix_sequence ) ;

    return CORK::matrix_valued_function::make_matrix_polynomial( polybasis, coefs );
  } // nonlinear_system()

} // namespace CORK

#endif

#endif
