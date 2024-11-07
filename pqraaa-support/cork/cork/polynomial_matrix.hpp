//  (C) Copyright Karl Meerbergen & Dries De Samblanx, 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_polynomial_matrix_hpp
#define cork_polynomial_matrix_hpp

#include <cork/linear_solver/user_defined.hpp>
#include <cork/utility/ref.hpp>
#include <cork/coefficient_matrices/matrices_by_functions.hpp>
#include <cork/matrix_valued_function/nonlinear_matrix.hpp>
#include <cork/basis/monomial.hpp>
#include <cork/concept/arithmetic.hpp>
#include <cork/concept/multiply_add.hpp>
#include <cork/concept/solve.hpp>

namespace CORK {

#ifdef CORK_USE_CONCEPTS
  template <typename ValueType, typename MultiplyAdd, typename Solve>
        requires CORK::Arithmetic<ValueType> && CORK::Solve< typename deref_type<Solve>::type,ValueType> && CORK::MultiplyAdd<typename deref_type<MultiplyAdd>::type,ValueType>
#else
  template <typename ValueType, typename MultiplyAdd, typename Solve>
#endif
  auto polynomial_matrix ( int n, int degree, MultiplyAdd multiply_add, Solve solve ) {
    assert( n>=1 ) ;
    assert( degree>=1 ) ;
    CORK::basis::monomial<>                                                      basis( degree );
    coefficient_matrices::matrices_by_functions< ValueType, MultiplyAdd, Solve > coefs( n, degree, multiply_add, solve ) ;

    return CORK::matrix_valued_function::nonlinear_matrix( basis, coefs );
  }


#ifdef CORK_USE_CONCEPTS
  template <typename ValueType, typename MultiplyAdd, typename Solve, typename Poly>
        requires CORK::Arithmetic<ValueType> && CORK::Solve<Solve,ValueType> && CORK::MultiplyAdd<MultiplyAdd,ValueType>
#else
  template <typename ValueType, typename MultiplyAdd, typename Solve>
#endif
  auto polynomial_matrix ( int n, Poly const& poly, MultiplyAdd const& multiply_add, Solve& solve ) {
    assert( n>=1 ) ;
    static_assert( !std::is_integral< Poly >::value, "CORK::polynomial_matrix: Poly cannot be integral") ;
    coefficient_matrices::matrices_by_functions< ValueType, MultiplyAdd, Solve > coefs( n, poly.num_terms()-1, multiply_add, solve ) ;

    return CORK::matrix_valued_function::nonlinear_matrix( poly, coefs );
  }
} // namespace CORK

// Insert for dense or sparse matrices if included
#if defined(cork_dense_hpp) || defined(cork_sparse_hpp)

namespace CORK {

  template <typename MatrixSequence>
  auto polynomial_matrix ( MatrixSequence const& matrices ) {
    CORK::basis::monomial<>                                   basis( matrices.size()-1 );
    coefficient_matrices::any_glas< MatrixSequence const& >   coefs( matrices ) ;

    return CORK::matrix_valued_function::nonlinear_matrix( basis, coefs );
  }
} // namespace CORK

#endif

#endif
