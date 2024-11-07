//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_sv_raaa_hpp
#define cork_matrix_valued_function_sv_raaa_hpp

#include <cork/matrix_valued_function/matrix_polynomial.hpp>
#include <cork/matrix_valued_function/sv_aaa.hpp>
//#include <cork/matrix_valued_function/matrix_polynomial_with_shadow.hpp>
#include <cork/matrix_valued_function/difference.hpp>
#include <cork/exception/option_out_of_bounds.hpp>
#include <cork/approximation/sv_aaa_real.hpp>
#include <cork/approximation/sv_aaa_poles_real.hpp>
#include <cork/basis/union_of_bases.hpp>
#include <cork/basis/barycentric_rational_real_strong.hpp>
#include <cork/options/check_poles.hpp>
#include <cork/options/compute_poles.hpp>
#include <cork/options/test_approximation_error.hpp>
#include <cork/options/value_of.hpp>
#include <cork/coefficient_matrices/combined.hpp>
#include <type_traits>
#include <cmath>
#include <cassert>

namespace CORK { namespace matrix_valued_function {

  namespace sv_raaa_detail {
    template <typename NEP, typename Domain, typename NormEst, typename Options>
    decltype (auto) sv_raaa_approximate ( NEP const& nep, Domain const& domain, NormEst& norm_est, Options const& options ) {
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

      typedef typename CORK::options::aaa_stop_criterion<real_type>::enum_type                aaa_stop_enum;
      typedef CORK::options::aaa_algorithm::enum_type                                         aaa_alg_enum;
      aaa_alg_enum aaa_alg = options::value_of<options::aaa_algorithm>(options);
      aaa_stop_enum aaa_stop = options::value_of<options::aaa_stop_criterion<real_type>>(options).choice();
      if (aaa_stop!=aaa_stop_enum::MAX) norm_est.activate() ;
      glas2::range r_aaa( nep.basis().basis_1().num_terms(), nep.basis().basis_1().num_terms()+nep.basis().basis_2().num_terms() ) ;

      // Scale function norm if needed for stop criterion
      sv_aaa_err_est<value_type,NormEst> err_est( norm_est, r_aaa, fun_sequence.num_terms(), aaa_stop ) ;
      err_est.scale_functions( test_set, test_values ) ;

      if (aaa_alg==aaa_alg_enum::SVD_SV_AAA) {
        std::cout << "CORK: AAA (real): SVD_SV_AAA is not yet implemented for real functions. Using SV_AAA instead." << std::endl ;
      }
      auto repr = approximation::SV_AAA_real( test_set, test_values, err_est, options ) ;

      err_est.unscale_functions( repr.coefficients() ) ;

      return repr ;
    }
  } // namespace sv_raaa_detail


  //
  // NonlinearMatrix must have the form MatrixPolynomial< basis::union_of_functions<Basis, basis::functions<T> >, CoefficientMatrices >
  // e.g., created with make_nonlinear_matrix().
  //
  template <typename NonlinearMatrix, typename AAANonlinearMatrix, typename Domain, typename NormEst, typename Options>
  decltype (auto) sv_raaa( NonlinearMatrix const& nep, AAANonlinearMatrix const& aaa_nep, Domain const& domain, NormEst& norm_est, Options const& options ) {
    //typedef T                            value_type ;
    //typedef typename NonlinearMatrix::basis_type::template value_type_for< typename Domain::value_type> value_type ;
    typedef typename Domain::value_type           value_type ;
    typedef decltype( std::abs( value_type() ) )  real_value_type ;

    auto const& monomial_basis = nep.basis().basis_1() ;
    auto const& nep_nonlinear = nep.basis().basis_2() ;

    auto aaa_approx = sv_raaa_detail::sv_raaa_approximate( aaa_nep, domain, norm_est, options );

    // Compute poles
    if (options::value_of<options::compute_poles>(options)) {
      auto poles = approximation::sv_aaa_poles_real( aaa_approx ) ;
      std::cout << "Poles " << poles << std::endl ;
    }

    /*// Sort nodes and weights so that real nodes come first.
    glas2::shared_vector< value_type > nodes( aaa_approx.nodes().size() ) ;
    glas2::shared_vector< value_type > weights( aaa_approx.nodes().size() ) ;
    glas2::matrix< real_value_type, glas2::row_major > coefficients( aaa_approx.coefficients().num_rows(), aaa_approx.coefficients().num_columns() ) ;
    int i_real_end = 0;
    int i_complex_end = aaa_approx.nodes().size() ;
    for (int i=0; i<aaa_approx.nodes().size(); ++i) {
      if (aaa_approx.nodes()(i).imag()==0.0) {
        nodes( i_real_end ) = aaa_approx.nodes()(i) ;
        weights( i_real_end ) = aaa_approx.weights()(i) ;
        coefficients( i_real_end, glas2::all() ) = aaa_approx.coefficients()(i, glas2::all()) ;
        ++i_real_end ;
      } else {
        assert( i<aaa_approx.nodes().size()-1) ;
        assert( std::conj(aaa_approx.nodes()(i))==aaa_approx.nodes()(i+1) ) ;
        nodes( i_complex_end-2 ) = aaa_approx.nodes()(i) ;
        nodes( i_complex_end-1 ) = aaa_approx.nodes()(i+1) ;
        weights( i_complex_end-2 ) = aaa_approx.weights()(i) ;
        weights( i_complex_end-1 ) = aaa_approx.weights()(i+1) ;
        coefficients( i_complex_end-2, glas2::all() ) = aaa_approx.coefficients()(i, glas2::all()) ;
        coefficients( i_complex_end-1, glas2::all() ) = aaa_approx.coefficients()(i+1, glas2::all()) ;
        i_complex_end -= 2 ;
        ++i ;
      }
    }
    assert( i_real_end==i_complex_end ) ;*/
    auto nodes( aaa_approx.nodes() ) ;
    auto weights( aaa_approx.weights() ) ;
    auto coefficients( aaa_approx.coefficients() ) ;

    //std::cout << "nodes " << nodes << std::endl ;
    //std::cout << "weights " << weights << std::endl ;

    // Make state_space basis
    typedef glas2::shared_vector< value_type > vector_type ;
    vector_type nodes_c( copy(nodes) ) ;
    vector_type weights_c( copy(weights) ) ;
    basis::barycentric_rational_real_strong< vector_type, vector_type > state_space_basis( weights_c, nodes_c ) ;

    auto uni_basis = basis::make_union_of_bases( monomial_basis, state_space_basis ) ;

    // Make combinations of coefficient matrices
    // The first (1,1) block is the identity matrix corresponding to the polynomial part
    // The combinations in the (2,2) block are equal to the function values of AAA multiplied with the weights, for each nonlinear functions, i.e., for each row of 'combinations'
    typedef glas2::shared_matrix< real_value_type > real_matrix_type ;
    real_matrix_type combinations( nep.coefficient_matrices().num_matrices(), uni_basis.num_terms() ) ;
    fill( combinations, 0.0 ) ;

    glas2::range basis_range(0,monomial_basis.num_terms()) ;
    // (1,1) block
    combinations( basis_range, basis_range ) = glas2::identity_matrix< real_value_type >( basis_range.size(), basis_range.size() ) ;

    // (2,2) block: the weights are taken into the linearization, not the matrices.
    combinations( glas2::range(monomial_basis.num_terms(),monomial_basis.num_terms()+nep_nonlinear.num_terms()), glas2::range_from_end(monomial_basis.num_terms(),0) ) = transpose( coefficients ) ;

    // (3,3) block: the weights are taken into the linearization, not the matrices.
    //combinations( glas2::range_from_end(monomial_basis.num_terms()+nep_nonlinear.num_terms(),0), glas2::range_from_end(monomial_basis.num_terms()+nodes.size(),0) ) = transpose( coefficients ) ;

    // Make the coefficient matrices for the AAA representation.
    CORK::coefficient_matrices::combined< typename NonlinearMatrix::coefficient_matrices_type, real_matrix_type > combined( nep.coefficient_matrices(), combinations ) ;

    // Make Rational eigenvalue problem, using a union of bases
    matrix_polynomial< decltype(uni_basis), decltype(combined) > pep( uni_basis, combined ) ;

    // Make Rational eigenvalue problem, using a union of bases
    CORK::coefficient_matrices::combined< typename AAANonlinearMatrix::coefficient_matrices_type, real_matrix_type > aaa_combined( aaa_nep.coefficient_matrices(), combinations ) ;
    matrix_polynomial< decltype(uni_basis), decltype(aaa_combined) > aaa_pep( uni_basis, aaa_combined ) ;

    if (options::value_of<options::debug_level>(options)>2) {
      std::cout << "Support points " << nodes_c << std::endl ;
    }

    if (options::value_of<options::test_approximation_error>(options)) {
       std::cout << "Number of terms of the polynomial matrix: " << pep.basis().num_terms() << std::endl ;
       std::cout << "Error on the AAA approximation " << CORK::matrix_valued_function::difference( aaa_nep, aaa_pep, domain, options ) << std::endl ;
    }

    /*
    if (options::value_of<options::test_approximation_error>(options)) {
       std::cout << "Number of terms of the polynomial matrix: " << pep.basis().num_terms() << std::endl ;
       CORK::coefficient_matrices::combined< typename AAANonlinearMatrix::coefficient_matrices_type const&, real_matrix_type > combined( aaa_nep.coefficient_matrices(), combinations ) ;
       auto uni_basis = basis::make_union_of_bases( aaa_nep.basis().basis_1(), state_space_basis ) ;
       matrix_polynomial< decltype(uni_basis), decltype(combined) > aaa_pep( uni_basis, combined ) ;
       std::cout << "Error on the AAA approximation " << CORK::matrix_valued_function::difference( aaa_nep, aaa_pep, domain, options ) << std::endl ;
    }*/

    //return make_matrix_polynomial_with_shadow( nep, pep ) ;
    return pep ;
    //return pep ;
  } // sv_raaa()


  template <typename NonlinearMatrix, typename Domain, typename NormEst>
  decltype (auto) sv_raaa ( NonlinearMatrix const& nep, Domain const& domain, NormEst& norm_est ) {
    return sv_raaa( nep, domain, norm_est, std::tuple<>() ) ;
  } // sv_aaa()


} } // namespace CORK::matrix_valued_function

#endif
