//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_sv_raaa_lsn_hpp
#define cork_matrix_valued_function_sv_raaa_lsn_hpp

#include <cork/matrix_valued_function/sv_raaa_ls.hpp>
//#include <cork/matrix_valued_function/matrix_polynomial_with_shadow.hpp>
#include <cork/approximation/sv_rat_accumulated_real_least_squares.hpp>

namespace CORK { namespace matrix_valued_function {

  //
  // NonlinearMatrix must have the form MatrixPolynomial< basis::union_of_functions<Basis, basis::functions<T> >, CoefficientMatrices >
  // e.g., created with make_nonlinear_matrix().
  //
  template <typename NonlinearMatrix, typename AAANonlinearMatrix, typename Domain, typename NormEst, typename Options>
  decltype (auto) sv_raaa_lsn( NonlinearMatrix const& nep, AAANonlinearMatrix const& aaa_nep, Domain const& domain, NormEst& norm_est, Options const& options ) {
    //typedef T                            value_type ;
    //typedef typename NonlinearMatrix::basis_type::template value_type_for< typename Domain::value_type> value_type ;
    typedef typename Domain::value_type           value_type ;
    typedef decltype( std::abs( value_type() ) )  real_value_type ;

    auto const& monomial_basis = nep.basis().basis_1() ;
    auto const& nep_nonlinear = nep.basis().basis_2().naked_functions() ;

    // Make united basis for rational approximation
    auto ls_approx = sv_raaa_ls_detail::sv_raaa_approximate( approximation::SV_rat_accumulated_real_least_squares_selector(), aaa_nep, domain, norm_est, options );
    auto uni_basis = basis::make_union_of_bases( monomial_basis, ls_approx.basis() ) ;

    auto coefficients( ls_approx.coefficients() ) ;

    //std::cout << "coefficients " << coefficients << std::endl ;

    // Make combinations of coefficient matrices
    // The first (1,1) block is the identity matrix corresponding to the polynomial part
    // The combinations in the (2,2) block are equal to the function values of AAA multiplied with the weights, for each nonlinear functions, i.e., for each row of 'combinations'
    typedef glas2::shared_matrix< real_value_type > real_matrix_type ;
    real_matrix_type combinations( nep.coefficient_matrices().num_matrices(), uni_basis.num_terms() ) ;
    fill( combinations, 0.0 ) ;

    glas2::range basis_range(0,monomial_basis.num_terms()) ;
    // (1,1) block
    combinations( basis_range, basis_range ) = glas2::identity_matrix< real_value_type >( basis_range.size(), basis_range.size() ) ;

    // Terms A_i * f_i
    // (2,2) block: the weights are taken into the linearization, not the matrices.
    combinations( glas2::range(monomial_basis.num_terms(),monomial_basis.num_terms()+nep_nonlinear.num_terms()), 0 ) = coefficients(0,glas2::all()) ;
    combinations( glas2::range(monomial_basis.num_terms(),monomial_basis.num_terms()+nep_nonlinear.num_terms()), glas2::range(monomial_basis.num_terms(),uni_basis.num_terms()) ) = transpose( coefficients(glas2::range_from_end(1,0),glas2::all() ) ) ;

    // Terms z * f_i * B_i
    // (3,3) block: the weights are taken into the linearization, not the matrices.
    auto weights = ls_approx.basis().weights() ;
    auto poles = ls_approx.basis().poles() ;
    static_assert( std::is_arithmetic< typename decltype(weights)::value_type >::value ) ;

    combinations( glas2::range_from_end(monomial_basis.num_terms()+nep_nonlinear.num_terms(),0), 1 ) = coefficients(0,glas2::all()) ;

    for (int i=0; i<poles.size(); ++i) {
      if (poles(i).imag()!=0.) {
        combinations( glas2::range_from_end(monomial_basis.num_terms()+nep_nonlinear.num_terms(),0), monomial_basis.num_terms()+i )
                = std::real(poles(i))*coefficients(i+1,glas2::all()) - glas2::imag(poles(i))*coefficients(i+2,glas2::all()) ;
        combinations( glas2::range_from_end(monomial_basis.num_terms()+nep_nonlinear.num_terms(),0), monomial_basis.num_terms()+i+1 )
                = glas2::real(poles(i))*coefficients(i+2,glas2::all()) + glas2::imag(poles(i))*coefficients(i+1,glas2::all()) ;
        ++i ;
      } else {
        combinations( glas2::range_from_end(monomial_basis.num_terms()+nep_nonlinear.num_terms(),0), monomial_basis.num_terms()+i )
                = glas2::real(poles(i))*coefficients(i+1,glas2::all()) ;
      }
    }

    typename decltype(weights)::value_type previous_weight = 1.0 ;
    int previous_1 = 0 ;
    for (int i=0; i<ls_approx.basis().num_real(); ++i) {
      combinations( glas2::range_from_end(monomial_basis.num_terms()+nep_nonlinear.num_terms(),0), previous_1 )
               += weights(i)/previous_weight*coefficients(i+1,glas2::all()) ;
      previous_1 = monomial_basis.num_terms()+i ;
      previous_weight = weights(i) ;
    }
    if (ls_approx.basis().num_real()<poles.size()) {
      combinations( glas2::range_from_end(monomial_basis.num_terms()+nep_nonlinear.num_terms(),0), previous_1 )
               += 2.*weights(ls_approx.basis().num_real())/previous_weight*coefficients(ls_approx.basis().num_real()+1,glas2::all()) ;

      previous_1 = monomial_basis.num_terms()+ls_approx.basis().num_real() ;
      previous_weight = weights(ls_approx.basis().num_real()) ;
    }

    for (int i=ls_approx.basis().num_real()+2; i<poles.size(); i+=2) {
      combinations( glas2::range_from_end(monomial_basis.num_terms()+nep_nonlinear.num_terms(),0), previous_1 )
              += weights(i)/previous_weight*coefficients(i+1,glas2::all()) ;
      combinations( glas2::range_from_end(monomial_basis.num_terms()+nep_nonlinear.num_terms(),0), previous_1+1 )
              += weights(i)/previous_weight*coefficients(i+2,glas2::all()) ;
      previous_1 = monomial_basis.num_terms()+i ;
      previous_weight = weights(i) ;
    }

    /*
    typename decltype(weights)::value_type previous_weight = 1.0 ;
    int previous_1 = 0 ;
    int previous_2 = 0 ;
    int last_pole = poles.size() ; if (poles(last_pole-1).imag()!=0) -- last_pole ;
    for (int i=0; i<last_pole; ++i) {
      if (poles(i).imag()!=0.) {
        combinations( glas2::range_from_end(monomial_basis.num_terms()+nep_nonlinear.num_terms(),0), previous_1 )
                += weights(i)/previous_weight*coefficients(i+1,glas2::all()) ;
        combinations( glas2::range_from_end(monomial_basis.num_terms()+nep_nonlinear.num_terms(),0), previous_2 )
                += weights(i)/previous_weight*coefficients(i+2,glas2::all()) ;
        previous_1 = monomial_basis.num_terms()+i ;
        previous_2 = monomial_basis.num_terms()+i+1 ;
      } else {
        combinations( glas2::range_from_end(monomial_basis.num_terms()+nep_nonlinear.num_terms(),0), previous_1 )
                += weights(i)/previous_weight*coefficients(i+1,glas2::all()) ;
        previous_1 = previous_2 = monomial_basis.num_terms()+i ;
      }
      previous_weight = weights(i) ;
    }

    /*
    int i1 = 1 ; if (poles(0).imag()!=0.) ++i1 ;
    combinations( glas2::range_from_end(monomial_basis.num_terms()+nep_nonlinear.num_terms(),0), 0 )
                += i1*weights(0)*coefficients(1,glas2::all()) ;

    for (int i=i1; i<poles.size(); ++i) {
      if (poles(i).imag()!=0.) {
        combinations( glas2::range_from_end(monomial_basis.num_terms()+nep_nonlinear.num_terms(),0), monomial_basis.num_terms()+i-i1 )
                += weights(i)/weights(i-i1)*coefficients(i+1,glas2::all()) ;
        combinations( glas2::range_from_end(monomial_basis.num_terms()+nep_nonlinear.num_terms(),0), monomial_basis.num_terms()+i-i1+1 )
                += weights(i)/weights(i-1)*coefficients(i+2,glas2::all()) ;
        ++i ;
      } else {
        combinations( glas2::range_from_end(monomial_basis.num_terms()+nep_nonlinear.num_terms(),0), monomial_basis.num_terms()+i-i1 )
                += weights(i)/weights(i-i1)*coefficients(i+1,glas2::all()) ;
      }
    }*/
    std::cout << combinations << std::endl ;

    // Make the coefficient matrices for the AAA representation.
    CORK::coefficient_matrices::combined< typename NonlinearMatrix::coefficient_matrices_type, real_matrix_type > combined( nep.coefficient_matrices(), combinations ) ;

    // Make Rational eigenvalue problem, using a union of bases
    matrix_polynomial< decltype(uni_basis), decltype(combined) > pep( uni_basis, combined ) ;

    // Make Rational eigenvalue problem, using a union of bases
    CORK::coefficient_matrices::combined< typename AAANonlinearMatrix::coefficient_matrices_type, real_matrix_type > aaa_combined( aaa_nep.coefficient_matrices(), combinations ) ;
    matrix_polynomial< decltype(uni_basis), decltype(aaa_combined) > aaa_pep( uni_basis, aaa_combined ) ;

    if (options::value_of<options::debug_level>(options)>2) {
      std::cout << "Poles " << ls_approx.basis().poles() << std::endl ;
      std::cout << "Weights " << ls_approx.basis().weights() << std::endl ;
    }

    if (options::value_of<options::test_approximation_error>(options)) {
       std::cout << "Number of terms of the polynomial matrix: " << pep.basis().num_terms() << std::endl ;
       std::cout << "Error on the AAA-LS approximation " << CORK::matrix_valued_function::difference( aaa_nep, aaa_pep, domain, options ) << std::endl ;
    }

    return pep ;
    //return pep ;
  } // sv_raaa_lsn()


  template <typename NonlinearMatrix, typename Domain, typename NormEst>
  decltype (auto) sv_raaa_lsn ( NonlinearMatrix const& nep, Domain const& domain, NormEst& norm_est ) {
    return sv_raaa_lsn( nep, domain, norm_est, std::tuple<>() ) ;
  } // sv_aaa_lsn()


} } // namespace CORK::matrix_valued_function

#endif
