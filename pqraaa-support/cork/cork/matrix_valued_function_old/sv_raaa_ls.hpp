//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_sv_raaa_ls_hpp
#define cork_matrix_valued_function_sv_raaa_ls_hpp

#include <cork/matrix_valued_function/matrix_polynomial.hpp>
#include <cork/matrix_valued_function/sv_aaa.hpp>
//#include <cork/matrix_valued_function/matrix_polynomial_with_shadow.hpp>
#include <cork/matrix_valued_function/difference.hpp>
#include <cork/exception/option_out_of_bounds.hpp>
#include <cork/approximation/sv_aaa_ls_real.hpp>
#include <cork/approximation/sv_aaa_filtered_ls_real.hpp>
#include <cork/approximation/sv_rat_real_least_squares.hpp>
#include <cork/basis/union_of_bases.hpp>
#include <cork/basis/partial_fractions_real.hpp>
#include <cork/options/sample_functions.hpp>
#include <cork/options/has_option.hpp>
#include <cork/options/check_poles.hpp>
#include <cork/options/compute_poles.hpp>
#include <cork/options/filter_poles.hpp>
#include <cork/options/test_approximation_error.hpp>
#include <cork/options/value_of.hpp>
#include <cork/coefficient_matrices/combined.hpp>
#include <type_traits>
#include <cmath>
#include <cassert>

namespace CORK { namespace matrix_valued_function {

  namespace sv_raaa_ls_detail {

    template <bool Separate>
    struct compute_samples {
    } ;

    template <>
    struct compute_samples<true> {
      template <typename TestSet, typename TestValues, typename Options>
      static auto apply( TestSet const& test_set, TestValues const& test_values, Options const& options ) {
        auto aaa_fun_sequence = options::value_of< options::sample_functions<typename TestValues::value_type> >( options ) ;
        glas2::shared_matrix<typename TestValues::value_type> aaa_values( test_set.size(), aaa_fun_sequence.num_terms() ) ;
        glas2::vector<typename TestValues::value_type> temp( aaa_values.num_columns() ) ;

        for (typename TestSet::size_type i=0; i<test_set.size(); ++i) {
          aaa_fun_sequence.evaluate( test_set(i), temp ) ;
          aaa_values(i,glas2::all()) = temp ;
          for (int j=0; j<aaa_values.num_columns(); ++j) {
            if (std::isnan(std::abs(test_values(i,j))) || std::isinf(std::abs(test_values(i,j)))) {
              if (options::value_of<options::debug_level>(options)>1) {
                std::cout << "Function values in the discrete points\n" ;
                std::cout << "Argument ->  Function values:\n" ;
                std::cout << test_set(i) << " -> " << test_values(i, j) << std::endl ;
              }
              throw exception::rational_approximation("Functions evaluate to Inf or NaN") ;
            }
          }
        }

        return aaa_values ;
      }
    } ;

    template <>
    struct compute_samples<false> {
      template <typename TestSet, typename TestValues, typename Options>
      static auto apply( TestSet const& test_set, TestValues const& test_values, Options const& options ) {
        return test_values ;
      }
    } ;

    template <typename AAA_LS_Selector, typename NEP, typename Domain, typename NormEst, typename Options>
    decltype (auto) sv_raaa_approximate ( AAA_LS_Selector const& AAA_LS_selector, NEP const& nep, Domain const& domain, NormEst& norm_est, Options const& options ) {
      // Construct AAA approximation
      // The functions are in the second part of the union_of_functions.

      typedef typename NEP::basis_type::template value_type_for< typename Domain::value_type> value_type ;
      typedef decltype(std::abs(value_type()))                                                real_type ;

      // Discretize the domain
      // Construct the test set
      auto test_set_ini = domain.discretize( std::max( 2*options::value_of<options::max_degree>(options), options::value_of<options::number_of_sample_points>(options)) ) ;
      auto test_set = approximation::unique_test_set( test_set_ini ) ;

      // Evaluate and scale the functions
      auto const& fun_sequence = nep.basis().basis_2() ;
      glas2::shared_matrix<value_type> test_values( test_set.size(), fun_sequence.num_terms() ) ;

      for (typename decltype(test_set)::size_type i=0; i<test_set.size(); ++i) {
        fun_sequence.evaluate( test_set(i), test_values(i,glas2::all()) ) ;
        for (int j=0; j<test_values.num_columns(); ++j) {
          if (std::isnan(std::abs(test_values(i,j))) || std::isinf(std::abs(test_values(i,j)))) {
            if (options::value_of<options::debug_level>(options)>1) {
              std::cout << "Function values in the discrete points\n" ;
              std::cout << "Argument ->  Function values:\n" ;
              std::cout << test_set(i) << " -> " << test_values(i, j) << std::endl ;
            }
            throw exception::rational_approximation("Functions evaluate to Inf or NaN") ;
          }
        }
      }

      // Get test values for AAA
      auto aaa_values = compute_samples< options::has_option< options::sample_functions<value_type>, Options >::value >::apply( test_set, test_values, options ) ;

      typedef typename CORK::options::aaa_stop_criterion<real_type>::enum_type                aaa_stop_enum;
      typedef CORK::options::aaa_algorithm::enum_type                                         aaa_alg_enum;
      aaa_alg_enum aaa_alg = options::value_of<options::aaa_algorithm>(options);
      aaa_stop_enum aaa_stop = options::value_of<options::aaa_stop_criterion<real_type>>(options).choice();
      if (aaa_stop!=aaa_stop_enum::MAX) norm_est.activate() ;
      //glas2::range r_aaa( nep.basis().basis_1().num_terms(), nep.basis().basis_1().num_terms()+nep.basis().basis_2().num_terms() ) ;
      // Assume that the last matrices correspond to the nonlinear terms
      glas2::range r_aaa( norm_est.size()-nep.basis().basis_2().num_terms(), norm_est.size() ) ;

      // Scale function norm if needed for stop criterion
      sv_aaa_err_est<value_type,NormEst> err_est( norm_est, r_aaa, fun_sequence.num_terms(), aaa_stop ) ;
      err_est.scale_functions( test_set, test_values ) ;

      if (aaa_alg==aaa_alg_enum::SVD_SV_AAA) {
        std::cout << "CORK: AAA (real): SVD_SV_AAA is not yet implemented for real functions. Using SV_AAA instead." << std::endl ;
      }
      auto filter = [&options]( CORK::vector< value_type > v ) {
        int n = options::value_of< options::filter_poles<value_type> >(options)( v ) ;
        return v( glas2::range(0,n) ) ;
      } ;
      auto repr = AAA_LS_selector( test_set, aaa_values, test_values, err_est, filter, nep.basis().basis_1().grade(), options ) ;

      err_est.unscale_functions( repr.coefficients() ) ;

      return repr ;
    }
  } // namespace sv_raaa_ls_detail


  //
  // NonlinearMatrix must have the form MatrixPolynomial< basis::union_of_functions<Basis, basis::functions<T> >, CoefficientMatrices >
  // e.g., created with make_nonlinear_matrix().
  //
  template <typename AAA_LS_Selector, typename NonlinearMatrix, typename AAANonlinearMatrix, typename Domain, typename NormEst, typename Options>
  decltype (auto) sv_raaa_ls( AAA_LS_Selector const& AAA_LS_selector, NonlinearMatrix const& nep, AAANonlinearMatrix const& aaa_nep, Domain const& domain, NormEst& norm_est, Options const& options ) {
    //typedef T                            value_type ;
    //typedef typename NonlinearMatrix::basis_type::template value_type_for< typename Domain::value_type> value_type ;
    typedef typename Domain::value_type           value_type ;
    typedef decltype( std::abs( value_type() ) )  real_value_type ;

    auto const& monomial_basis = nep.basis().basis_1() ;
    auto const& nep_nonlinear = nep.basis().basis_2() ;

    // Make united basis for rational approximation
    auto ls_approx = sv_raaa_ls_detail::sv_raaa_approximate( AAA_LS_selector, aaa_nep, domain, norm_est, options );
    //auto uni_basis = basis::make_union_of_bases( monomial_basis, ls_approx.basis() ) ;
    auto uni_basis = ls_approx.basis()  ;

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
    combinations( glas2::range(monomial_basis.num_terms(),monomial_basis.num_terms()+nep_nonlinear.num_terms()), glas2::all() ) = transpose( coefficients ) ;
    std::cout << combinations << std::endl ;

    // Make the coefficient matrices for the AAA representation.
    CORK::coefficient_matrices::combined< typename NonlinearMatrix::coefficient_matrices_type, real_matrix_type > combined( nep.coefficient_matrices(), combinations ) ;

    // Make Rational eigenvalue problem, using a union of bases
    matrix_polynomial< decltype(uni_basis), decltype(combined) > pep( uni_basis, combined ) ;

    // Make Rational eigenvalue problem, using a union of bases
    CORK::coefficient_matrices::combined< typename AAANonlinearMatrix::coefficient_matrices_type, real_matrix_type > aaa_combined( aaa_nep.coefficient_matrices(), combinations ) ;
    matrix_polynomial< decltype(uni_basis), decltype(aaa_combined) > aaa_pep( uni_basis, aaa_combined ) ;

    if (options::value_of<options::test_approximation_error>(options)) {
       std::cout << "Number of terms of the polynomial matrix: " << pep.basis().num_terms() << std::endl ;
       std::cout << "Error on the AAA-LS approximation " << CORK::matrix_valued_function::difference( aaa_nep, aaa_pep, domain, options ) << std::endl ;
    }

    return pep ;
    //return pep ;
  } // sv_raaa_ls()

  // Extended AAA
  template <typename NonlinearMatrix, typename AAANonlinearMatrix, typename Domain, typename NormEst, typename Options>
  decltype (auto) sv_raaa_ls( NonlinearMatrix const& nep, AAANonlinearMatrix const& aaa_nep, Domain const& domain, NormEst& norm_est, Options const& options ) {
    return sv_raaa_ls( approximation::SV_AAA_LS_real_selector(), nep, aaa_nep, domain, norm_est, options ) ;
  } // sv_aaa_ls()

  template <typename NonlinearMatrix, typename Domain, typename NormEst>
  decltype (auto) sv_raaa_ls ( NonlinearMatrix const& nep, Domain const& domain, NormEst& norm_est ) {
    return sv_raaa_ls( approximation::SV_AAA_filtered_LS_real_selector(), nep, nep, domain, norm_est, std::tuple<>() ) ;
  } // sv_aaa_ls()


  // Filtered Extended AAA
  template <typename NonlinearMatrix, typename AAANonlinearMatrix, typename Domain, typename NormEst, typename Options>
  decltype (auto) sv_raaa_filtered_ls( NonlinearMatrix const& nep, AAANonlinearMatrix const& aaa_nep, Domain const& domain, NormEst& norm_est, Options const& options ) {
    return sv_raaa_ls( approximation::SV_AAA_filtered_LS_real_selector(), nep, aaa_nep, domain, norm_est, options ) ;
  } // sv_aaa_ls()



} } // namespace CORK::matrix_valued_function

#endif
