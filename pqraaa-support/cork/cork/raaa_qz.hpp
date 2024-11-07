//  (C) Copyright Karl Meerbergen 2023.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_raaa_qz_hpp
#define cork_raaa_qz_hpp

#include <cork/eigs/info.hpp>
#include <cork/eigs/qz_pair.hpp>
#include <cork/eigs/do_qz.hpp>
#include <cork/linearization/union_of_bases.hpp>
#include <cork/linearization/barycentric_rational_real_strong.hpp>
#include <cork/linearization/combined.hpp>
#include <cork/linearization/cork_linearization.hpp>
#include <cork/linearization/monomial.hpp>
#include <cork/matrix_valued_function/select_nep_from_pair.hpp>
#include <cork/matrix_valued_function/nonlinear_matrix_with_shadow.hpp>
#include <cork/matrix_valued_function/linear_solver.hpp> // For dense and sparse solvers
#include <cork/matrix_valued_function/user_defined.hpp> // For user defined solve
#include <cork/matrix_valued_function/sv_raaa.hpp>
#include <cork/matrix_valued_function/norm_est.hpp>
#include <cork/matrix_valued_function/difference.hpp>
#include <cork/options/value_of.hpp>
#include <cork/eigs/stop_criterion.hpp>
#include <cork/options/debug_level.hpp>
#include <cork/options/shift_selection.hpp>
#include <cork/options/handle_as_real_problem.hpp>
#include <cork/options/has_option.hpp>
#include <cork/options/use_pep_linear_solver.hpp>
#include <cork/options/test_approximation_error.hpp>
#include <cork/utility/select_if.hpp>
#include <glas2/bindings/vector.hpp>
#include <tuple>
#include <typeinfo>

namespace CORK {

 /**
  * The actual aaa_cork algorithm.
  * 
  * Returns amongst other things the eigenvalues. 
  * 
  * @param NEP
  * @param aaaDomain
  * @param eigenvalue_selector
  * @param n_wanted
  * @param Options
  */
  template <typename NEP, typename Domain, typename EigenvalueSelector, typename Options>
  decltype (auto) raaa_qz( NEP const& problem, Domain const& aaa_domain, EigenvalueSelector const& eig_selector, Options const& options ) {
    std::vector< std::string > warnings ;

    auto const& nep_krylov = matrix_valued_function::select_nep_from_pair< typename Domain::value_type, std::false_type >( CORK::deref(problem) ) ;
    typedef typename std::decay<decltype(nep_krylov)>::type nep_krylov_type ;

    typedef typename nep_krylov_type::template value_type_for< typename Domain::value_type >                                  value_type ;
    typedef decltype(std::abs(value_type()))                                                                                                    real_type ;
    typedef std::complex< real_type >                                                                                                           eig_value_type ;

    // norm estimation of matrix for weighted AAA approximation
    auto const& nep_eig = matrix_valued_function::select_nep_from_pair< eig_value_type, std::true_type >( CORK::deref(problem) ) ;
    matrix_valued_function::norm_est norm_est(nep_eig) ;

    auto const& nep_aaa = matrix_valued_function::select_nep_from_pair< typename Domain::value_type, std::true_type >( CORK::deref(problem) ) ;
    //matrix_valued_function::nonlinear_matrix_with_shadow PEP( CORK::matrix_valued_function::sv_aaa( nep_aaa, aaa_domain, norm_est, options ), nep_krylov ) ;
    auto [PEP, aaa_info] = CORK::matrix_valued_function::sv_raaa( nep_aaa, nep_aaa, aaa_domain, norm_est, options ) ;

    if (options::value_of<options::test_approximation_error>(options)) {
      matrix_valued_function::difference ( nep_aaa, PEP, aaa_domain, options ) ;
    }

    if (options::value_of<options::debug_level>(options)>0)
      std::cout << "Degree of polynomial matrix " << PEP.num_terms() << std::endl ;

    eigs::info information ;
    information.approximation = aaa_info ;
    auto linearization = CORK::linearization::make_cork_linearization( &PEP, information ) ;
    auto fill_handle = linearization.fill_handle() ;

    eigs::qz_pair<value_type> qz( linearization.num_rows() ) ;
    fill_handle.A( qz.A_ ) ;
    fill_handle.B( qz.B_ ) ;

    //Compute Schur decomposition

    return std::tuple< eigs::do_qz, decltype(qz), decltype(linearization), decltype(eig_selector) const&, decltype(information) >( eigs::do_qz(), std::move(qz), linearization, eig_selector, information );
  } // aaa_qz()

} // namespace CORK

#endif
