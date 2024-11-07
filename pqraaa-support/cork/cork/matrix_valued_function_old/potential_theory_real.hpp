//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_potential_theory_real_hpp
#define cork_matrix_valued_function_potential_theory_real_hpp

#include <cassert>
#include <cork/matrix_valued_function/potential_theory_options.hpp>
#include <cork/matrix_valued_function/matrix_polynomial.hpp>
#include <cork/matrix_valued_function/difference.hpp>
#include <cork/approximation/potential_theory.hpp>
#include <cork/approximation/potential_theory_real.hpp>
#include <cork/basis/union_of_bases.hpp>
#include <cork/coefficient_matrices/combined.hpp>
#include <type_traits>
#include <stdexcept>
#include <vector>


namespace CORK { namespace matrix_valued_function {


  //
  // NonlinearMatrix must have the form MatrixPolynomial< basis::union_of_functions<Basis, basis::functions<T> >, CoefficientMatrices >
  // e.g., created with make_nonlinear_matrix().
  //
  template <typename T, typename NonlinearMatrix, typename Sigma, typename Xi>
  decltype (auto) potential_theory_real ( NonlinearMatrix const& nep, Sigma const& sigma_domain, Xi const& xi_domain ) {
    return potential_theory_real<T>( nep, sigma_domain, xi_domain, matrix_valued_function::potential_theory_options< decltype(std::abs(T())) >() ) ;
  } // potential_theory()


  template <typename T, typename NonlinearMatrix, typename Sigma, typename Xi, typename Options>
  decltype (auto) potential_theory_real ( NonlinearMatrix const& nep, Sigma const& sigma_domain, Xi const& xi_domain, Options const& options ) {
    typedef T value_type ;
    typedef glas2::shared_vector< value_type > vector_type ;

    // Construct Rational Newton approximation
    vector_type sigma_points = sigma_domain.border().discretize( options.n_points ) ;
    vector_type xi_points    = xi_domain.border().discretize( options.n_points ) ;
    // Check that xi and sigma have no points in common.
    // TODO
    
    // Construct the Newton basis
    vector_type nodes( options.max_degree+1 ) ;
    vector_type poles( options.max_degree ) ;
    vector_type scaling( options.max_degree+1 ) ;
    CORK::basis::rational_newton< vector_type& > newton_basis( nodes, poles, scaling ) ;

    // The functions are in the second part of the union_of_functions.
    approximation::potential_theory_leja_bagby_real( sigma_points, xi_points, newton_basis ) ;
    approximation::potential_theory_scaling( sigma_points, newton_basis ) ;

    auto repr = approximation::potential_theory<value_type>( nep.basis().basis_2(), newton_basis, sigma_points, options ) ;

    // Make combinations of coefficient matrices
    // The first (1,1) block is the identity matrix corresponding to the polynomial part
    // The combinations in the (2,2) block are equal to the function values of AAA multiplied with the weights, for each nonlinear functions, i.e., for each row of 'combinations'
    typedef glas2::shared_matrix< value_type > matrix_type ;
    matrix_type combinations( nep.coefficient_matrices().num_matrices(), nep.basis().basis_1().num_terms()+repr.basis().num_terms()-1 ) ;
    fill( combinations, 0.0 ) ;

    glas2::range basis_range(0,nep.basis().basis_1().num_terms()) ;
    // (1,1) block
    combinations( basis_range, basis_range ) = glas2::identity_matrix< value_type >( basis_range.size(), basis_range.size() ) ;

    // (2,2) block: the weights are taken into the linearization, not the matrices.
    combinations( glas2::range_from_end(nep.basis().basis_1().num_terms(),0), 0 ) = repr.basis().scaling()(0) * repr.coefficients()(0,glas2::all()) ;
    combinations( glas2::range_from_end(nep.basis().basis_1().num_terms(),0), glas2::range_from_end(nep.basis().basis_1().num_terms(),0) ) = repr.basis().scaling()(0) * transpose( repr.coefficients()(glas2::range_from_end(1,0), glas2::all() ) ) ;
/*    for (int i=nep.basis().num_terms(); i<combinations.num_rows(); ++i) {
      combinations( i, glas2::range_from_end(nep.basis().num_terms(),0) ) = combinations( i, glas2::range_from_end(nep.basis().num_terms(),0) ) * potential_theory_approx.weights() ;
    }*/

    // Make the coefficient matrices for the NLEIGS representation.
    CORK::coefficient_matrices::combined< typename NonlinearMatrix::coefficient_matrices_type const&, matrix_type > combined( nep.coefficient_matrices(), combinations ) ;

    // Make Rational eigenvalue problem, using a union of bases
    auto uni_basis = basis::make_union_of_bases( nep.basis().basis_1(), repr.basis() ) ;
    matrix_polynomial< decltype(uni_basis), decltype(combined) > pep( uni_basis, combined ) ;

    if (options.test_error) {
       std::cout << "Degree of the matrix polynomial " << pep.basis().num_terms() << std::endl ;
       std::cout << "Error on the NLEIGS approximation " << CORK::matrix_valued_function::difference( pep, nep, sigma_domain ) << std::endl ;
    }

    return pep ;
  } // potential_theory()


} } // namespace CORK::matrix_valued_function


#endif
