//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_sv_raaa_strong_hpp
#define cork_matrix_valued_function_sv_raaa_strong_hpp

#include <cork/matrix_valued_function/matrix_polynomial.hpp>
#include <cork/matrix_valued_function/sv_raaa.hpp>
//#include <cork/matrix_valued_function/matrix_polynomial_with_shadow.hpp>
#include <cork/basis4cork/monomial_plus_barycentric_rational.hpp>
#include <cork/matrix_iterator/monomial_plus_barycentric_rational.hpp>
#include <cork/basis4cork/barycentric_rational_real.hpp>

namespace CORK { namespace matrix_valued_function {

  //
  // NonlinearMatrix must have the form MatrixPolynomial< basis::union_of_functions<Basis, basis::functions<T> >, CoefficientMatrices >
  // e.g., created with make_nonlinear_matrix().
  //
  template <typename NonlinearMatrix, typename AAANonlinearMatrix, typename Domain, typename NormEst, typename Options>
  decltype (auto) sv_raaa_strong( NonlinearMatrix const& nep, AAANonlinearMatrix const& aaa_nep, Domain const& domain, NormEst& norm_est, Options const& options ) {
    auto sv_raa_pep = sv_raaa( nep, aaa_nep, domain, norm_est, options ) ;

    // Make Rational eigenvalue problem, using a union of bases
    basis::monomial_plus_barycentric_rational new_basis( sv_raa_pep.basis().basis_1().grade(), sv_raa_pep.basis().basis_2() ) ;

    matrix_polynomial< decltype(new_basis), typename std::decay<decltype(sv_raa_pep.coefficient_matrices())>::type> pep( new_basis, sv_raa_pep.coefficient_matrices() ) ;

    return pep ;
    //return pep ;
  } // sv_raaa_strong()


} } // namespace CORK::matrix_valued_function

#endif
