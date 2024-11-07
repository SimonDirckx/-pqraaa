//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_sv_aaa_hpp
#define cork_matrix_valued_function_sv_aaa_hpp

#include <cassert>
#include <cork/matrix_valued_function/nonlinear_matrix.hpp>
#include <cork/matrix_valued_function/nonlinear_matrix_with_shadow.hpp>
#include <cork/matrix_valued_function/difference.hpp>
#include <cork/matrix_valued_function/sv_aaa_err_est.hpp>
#include <cork/approximation/fun2table.hpp>
#include <cork/approximation/sv_aaa_approximation.hpp>
#include <cork/approximation/svd_aaa.hpp>
#include <cork/approximation/svd_aaa_ls.hpp>
#include <cork/approximation/sv_aaa_max.hpp>
#include <cork/approximation/unique_test_set.hpp>
#include <cork/approximation/sv_aaa_2_partial_fractions.hpp>
#include <cork/options/has_option.hpp>
#include <cork/options/check_poles.hpp>
#include <cork/options/test_approximation_error.hpp>
#include <cork/options/aaa_algorithm.hpp>
#include <cork/options/aaa_stop_criterion.hpp>
#include <cork/options/value_of.hpp>
#include <cork/options/domain_distance_tolerance.hpp>
#include <cork/options/use_shadowed_nonlinear_matrix.hpp>
#include <cork/basis/union_of_bases.hpp>
#include <cork/basis/barycentric_rational_strong.hpp>
#include <cork/coefficient_matrices/combined.hpp>
#include <cork/exception/option_out_of_bounds.hpp>
#include <type_traits>
#include <stdexcept>
#include <cmath>
#include <vector>

namespace CORK { namespace matrix_valued_function {

  namespace sv_aaa_detail {

    template <typename PEP, typename NEP>
    auto shadow_or_pep( std::true_type, PEP const& pep, NEP const& nep ) {
      return nonlinear_matrix_with_shadow<PEP,NEP>( pep, nep ) ;
    }

    template <typename PEP, typename NEP>
    auto shadow_or_pep( std::false_type, PEP const& pep, NEP const& nep ) {
      return pep ;
    }

    template <typename AAAAlg, typename TestSet, typename FunVals, typename ErrEst, typename Options>
    decltype (auto) sv_or_svd_aaa( AAAAlg aaa_alg, TestSet const& test_set, FunVals const& fun_vals, ErrEst& err_est, Options const& options ) {
      switch(aaa_alg) {
        case AAAAlg::SV_AAA :
          return approximation::SV_AAA( test_set, fun_vals, err_est, options, false ) ;
        case AAAAlg::SVD_SV_AAA :
          return approximation::SVD_AAA( test_set, fun_vals, err_est, options, false ) ;
        default:
          throw exception::option_out_of_bounds("aaa_algorithm", "allowed values: SV_AAA | SVD_AAA");
          break;
      }
    } // sv_or_svd_aaa

    template <typename NEP, typename Domain, typename NormEst, typename Options>
    decltype (auto) sv_aaa_approximate ( NEP const& nep, Domain const& domain, NormEst& norm_est, Options const& options ) {
      // Construct AAA approximation
      // The functions are in the second part of the union_of_functions.

      typedef typename NEP::basis_type::template value_type_for< typename Domain::value_type> value_type ;
      typedef decltype(std::abs(value_type()))                                                real_type ;

      // Discretize the functions
      auto const& fun_sequence = nep.basis().basis_2() ;
      auto [test_set, test_values ] = approximation::fun2table( fun_sequence, domain, options ) ;

      typedef typename CORK::options::aaa_stop_criterion<real_type>::enum_type                aaa_stop_enum;
      typedef CORK::options::aaa_algorithm::enum_type                                         aaa_alg_enum;
      aaa_alg_enum aaa_alg = options::value_of<options::aaa_algorithm>(options);
      aaa_stop_enum aaa_stop = options::value_of<options::aaa_stop_criterion<real_type>>(options).choice();
      if (aaa_stop!=aaa_stop_enum::MAX) norm_est.activate() ;
      glas2::range r_aaa( nep.basis().basis_1().num_terms(), nep.basis().num_terms() ) ;

      // Scale function norm if needed for stop criterion
      sv_aaa_err_est<value_type,NormEst> err_est( norm_est, r_aaa, nep.basis().basis_2().num_terms(), aaa_stop ) ;
      err_est.scale_functions( test_set, test_values ) ;

      auto repr = sv_or_svd_aaa( aaa_alg, test_set, test_values, err_est, options ) ;
//      auto repr = AAA_LS_selector( test_set, aaa_values, test_values, err_est, filter, nep.basis().basis_1().grade(), options ) ;

      err_est.unscale_functions( repr.coefficients() ) ;

      return repr ;
    }
  } // namespace sv_aaa_detail


  //
  // NonlinearMatrix must have the form MatrixPolynomial< basis::union_of_functions<Basis, basis::functions<T> >, CoefficientMatrices >
  // e.g., created with make_nonlinear_matrix().
  //
  template <typename NonlinearMatrix, typename Domain, typename NormEst, typename Options>
  decltype (auto) sv_aaa( NonlinearMatrix const& nep, Domain const& domain, NormEst& err_est, Options const& options ) {
    typedef typename Domain::value_type                                          lambda_value_type ;
    typedef typename NonlinearMatrix::template value_type_for<lambda_value_type> value_type ;
    typedef decltype(std::abs(lambda_value_type()))                              real_type ;

    auto aaa_approx = sv_aaa_detail::sv_aaa_approximate( nep, domain, err_est, options );

    if (options::value_of<options::test_approximation_error>(options) && options::value_of<options::debug_level>(options)>1) {
      auto points = domain.discretize( options::value_of<options::number_of_sample_points>(options) ) ;
      glas2::vector<value_type> fun_values(aaa_approx.coefficients().num_columns()) ;
      for (int i=0; i<points.size(); ++i) {
        aaa_approx.eval(points(i), fun_values) ;
        std::cout << "AAA: " << points(i) << " --> " << fun_values << " " ;
        nep.basis().basis_2().evaluate( points(i), fun_values ) ;
        std::cout << " --> " << fun_values << std::endl ;
      }
    }

    if (aaa_approx.n()>0 && options::value_of<options::check_poles>(options)) {
      auto pf = approximation::sv_aaa_2_partial_fractions( aaa_approx ) ;
      if (options::value_of<options::debug_level>(options)>0) std::cout << "AAA poles " << pf.nodes() << std::endl ;
      int fail = 0 ;
      for (int i=0; i<pf.nodes().size(); ++i) {
        if (domain.distance( pf.nodes()(i) ) < options::value_of<options::domain_distance_tolerance<real_type>>(options)) {
          ++fail ;
          std::cerr << "Pole " << pf.nodes()(i) << " lies inside the domain" << std::endl ;
        }
      }
      if (fail>0) throw std::runtime_error( "There are AAA poles that lie inside the domain" ) ;
    }

    // Make state_space basis
    glas2::shared_vector<lambda_value_type> nodes( copy( aaa_approx.nodes() ) ) ;
    glas2::shared_vector<value_type> weights( copy( aaa_approx.weights() ) ) ;
    basis::barycentric_rational_strong< decltype(weights), decltype(nodes) > state_space_basis( weights, nodes ) ;

    // Make combinations of coefficient matrices
    // The first (1,1) block is the identity matrix corresponding to the polynomial part
    // The combinations in the (2,2) block are equal to the function values of AAA multiplied with the weights, for each nonlinear functions, i.e., for each row of 'combinations'
    typedef glas2::shared_matrix< value_type > matrix_type ;
    matrix_type combinations( nep.coefficient_matrices().num_matrices(), nep.basis().basis_1().num_terms()+nodes.size() ) ;
    fill( combinations, 0.0 ) ;

    glas2::range basis_range(0,nep.basis().basis_1().num_terms()) ;
    // (1,1) block
    combinations( basis_range, basis_range ) = glas2::identity_matrix< value_type >( basis_range.size(), basis_range.size() ) ;

    // (2,2) block: the weights are taken into the linearization, not the matrices.
    combinations( glas2::range_from_end(nep.basis().basis_1().num_terms(),0), glas2::range_from_end(nep.basis().basis_1().num_terms(),0) ) = transpose( aaa_approx.coefficients() ) ;
/*    for (int i=nep.basis().num_terms(); i<combinations.num_rows(); ++i) {
      combinations( i, glas2::range_from_end(nep.basis().num_terms(),0) ) = combinations( i, glas2::range_from_end(nep.basis().num_terms(),0) ) * aaa_approx.weights() ;
    }*/

    // Make the coefficient matrices for the AAA representation.
    CORK::coefficient_matrices::combined< typename NonlinearMatrix::coefficient_matrices_type const&, matrix_type > combined( nep.coefficient_matrices(), combinations ) ;

    // Make Rational eigenvalue problem, using a union of bases
    auto uni_basis = basis::make_union_of_bases( nep.basis().basis_1(), state_space_basis ) ;
    nonlinear_matrix< decltype(uni_basis), decltype(combined) > pep( uni_basis, combined ) ;

    if (options::value_of<options::test_approximation_error>(options)) {
       std::cout << "Degree of the matrix polynomial " << pep.basis().num_terms() << std::endl ;
       std::cout << "Error on the AAA approximation " << CORK::matrix_valued_function::difference( pep, nep, domain, options ) << std::endl ;
    }

    return std::tuple( sv_aaa_detail::shadow_or_pep( typename options::has_option<options::use_shadowed_nonlinear_matrix,Options>::type(), pep, nep )
                     , aaa_approx.info()
                     ) ;
    //return make_nonlinear_matrix_with_shadow( pep, nep ) ;
  } // sv_aaa()

  template <typename NonlinearMatrix, typename Domain, typename NormEst>
  decltype (auto) sv_aaa ( NonlinearMatrix const& nep, Domain const& domain, NormEst& norm_est ) {
    return sv_aaa( nep, domain, norm_est, std::tuple<>() ) ;
  } // sv_aaa()

} } // namespace CORK::matrix_valued_function

#endif
