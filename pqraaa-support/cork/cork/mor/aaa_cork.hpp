//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_aaa_cork_hpp
#define cork_aaa_cork_hpp

#ifndef CORK_CORK2
#include <cork/eigs/cork3.hpp>
#else
#include <cork/eigs/cork2.hpp>
#endif
#include <cork/eigs/cork_result.hpp>
#include <cork/eigs/info.hpp>
#include <cork/linearization/union_of_bases.hpp>
#include <cork/linearization/barycentric_rational.hpp>
#include <cork/linearization/combined.hpp>
#include <cork/linearization/cork_linearization.hpp>
#include <cork/linearization/monomial.hpp>
#include <cork/shift_generator/multiple_shifts.hpp>
#include <cork/matrix_valued_function/select_nep_from_pair.hpp>
#include <cork/matrix_valued_function/matrix_polynomial_with_shadow.hpp>
#include <cork/matrix_valued_function/sv_aaa.hpp>
#include <cork/matrix_valued_function/norm_est.hpp>
#include <cork/options/value_of.hpp>
#include <cork/eigs/stop_criterion.hpp>
#include <cork/options/debug_level.hpp>
#include <cork/options/initial_vector.hpp>
#include <cork/options/shift_selection.hpp>
#include <cork/utility/value_type_for.hpp>
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
  decltype (auto) aaa_cork( NEP const& problem, Domain const& aaa_domain, EigenvalueSelector const& eig_selector, Options const& options ) {
    std::vector< std::string > warnings ;

    auto shift_selector = produce_shift_selection(eig_selector,options);
    auto const& shifts = shift_selector.shifts() ;
    shift_generator::multiple_shifts< decltype(shifts) > shift_generator( shifts, shift_selector.multiplicity() ) ;

    auto const& nep_krylov = matrix_valued_function::select_nep_from_pair< typename decltype(shift_selector)::value_type, std::false_type >( problem ) ;
    typedef typename std::decay<decltype(nep_krylov)>::type nep_krylov_type ;

    typedef typename nep_krylov_type::template value_type_for< typename decltype(shift_selector)::value_type > value_type ;
    typedef decltype(std::abs(value_type()))                                                                   real_type ;
    typedef std::complex< real_type >                                                                          eig_value_type ;

    auto const& nep_eig = matrix_valued_function::select_nep_from_pair< eig_value_type, std::true_type >( problem ) ;
    typedef typename std::decay<decltype(nep_eig)>::type nep_eig_type ;
    matrix_valued_function::norm_est norm_est(nep_eig) ;

    auto PEP = make_matrix_polynomial_with_shadow( CORK::matrix_valued_function::sv_aaa( nep_krylov, aaa_domain, norm_est, options ), nep_krylov ) ;
    if (options::value_of<options::debug_level>(options)>0)
      std::cout << "Degree of rational AAA polynomial " << PEP.num_terms() << std::endl ;

    //matrix_valued_function::norm_est norm_est_pep(PEP) ;
    //eigs::stop_criterion<eig_value_type,decltype(PEP),options::stop_criterion<real_type>> stop_criterion( PEP, norm_est_pep, CORK::options::value_of<options::stop_criterion<real_type>>(options) ) ;
    eigs::stop_criterion<eig_value_type,nep_eig_type,options::stop_criterion<real_type>> stop_criterion( nep_eig, norm_est, CORK::options::value_of<options::stop_criterion<real_type>>(options) ) ;

    eigs::info information ;
    auto linearization = CORK::linearization::make_cork_linearization( PEP, information ) ;

    // Create the CORK triple.
    int n_wanted = eig_selector.n_wanted_max();
    CORK::krylov::cork_quadruple<value_type,decltype(options)> quadruple( linearization, n_wanted, options ) ;

    // Choose a random initial vector of full rank of the form
    // Q U^T where Q in R^{n\times 2} and U in R^2
    // via optie?
    // random 1 en functie?
    auto initial_vector = options::value_of<options::initial_vector<value_type>>(options);

    // Compute the Schur decomposition.
    CORK::eigs::cork( linearization, eig_selector, shift_generator, initial_vector, quadruple, stop_criterion, options );

    return CORK::eigs::cork_result(std::move(quadruple), linearization, eig_selector, information);
//    return CORK::eigs::make_cork_result( std::forward< decltype(quadruple) >(quadruple), 
//                                         std::forward< decltype(linearization) >(linearization),
//                                         std::forward< decltype(eig_selector) >(eig_selector), 
//                                         information )  ;

  } // aaa_cork()

} // namespace CORK

#endif
