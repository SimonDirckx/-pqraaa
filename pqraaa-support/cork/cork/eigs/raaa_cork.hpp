//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_eigs_raaa_cork_hpp
#define cork_eigs_raaa_cork_hpp

#ifdef cork_eigs_aaa_cork_hpp
#error raa_cork_hpp and aaa_cork.hpp cannot be included together
#endif

#include <cork/eigs/cork4.hpp>
#include <cork/eigs/info.hpp>
#include <cork/eigs/cork_eigen_pairs.hpp>
#include <cork/eigs/stop_criterion.hpp>
#include <cork/krylov/random_initial_vector.hpp>
#include <cork/linearization/union_of_bases.hpp>
#include <cork/linearization/cork_linearization.hpp>
#include <cork/linearization/barycentric_rational_real_strong.hpp>
#include <cork/linearization/monomial.hpp>
#include <cork/linearization/combined.hpp>
#include <cork/shift_generator/multiple_shifts.hpp>
#include <cork/matrix_valued_function/norm_est.hpp>
#include <cork/matrix_valued_function/sv_raaa.hpp>
#include <cork/eigenvalue_selector/inside_domain.hpp>
#include <cork/approximation/sv_aaa_real.hpp>
#include <cork/options/value_of.hpp>
#include <cork/options/shift_selection.hpp>
#include <cork/options/debug_level.hpp>
#include <cork/options/stop_criterion.hpp>
#include <glas2/bindings/vector.hpp>
#include <boost/numeric/bindings/std/vector.hpp>
#include <typeinfo>

namespace CORK { namespace eigs { namespace RAAA_CORK {

  template <typename T>
  struct eigen_system {
    eigen_system( int n, int n_cvgd )
    : eigen_values(n_cvgd)
    , eigen_vectors(n,n_cvgd)
    , shifts(0)
    {}

    glas2::vector<T> eigen_values ;
    glas2::matrix<T> eigen_vectors ;
    eigs::info       information ;
    glas2::vector<T> shifts ;
  } ;


  // Options should contain which AAA, which shift_generator to choose
/*  template <typename Scenario, typename T=typename Scenario::value_type>
  struct options {
    typedef decltype(std::abs(T())) real_type ;
    typedef T value_type ;

    options( Scenario const& shift_generator )
    : shift_generator( shift_generator )
    {}

    matrix_valued_function::aaa_options<real_type> aaa_options ;
    eigs::options<real_type>                       krylov_options ;
    Scenario                                       shift_generator ;
  } ;

  template <typename Scenario>
  decltype (auto) make_options( Scenario const& shift_generator ) { return options<Scenario>( shift_generator ) ; }
*/

  template <typename NEP, typename AAADomain, typename EigDomain, typename Options>
  decltype (auto) compute_eigenvalues_and_vectors( NEP const& problem, AAADomain const& aaa_domain, EigDomain const& eig_domain, int n_wanted, Options const& options ) {
    std::vector< std::string > warnings ;

    matrix_valued_function::norm_est norm_est(problem) ;

    auto PEP = CORK::matrix_valued_function::sv_raaa( problem, problem, aaa_domain, norm_est, options ) ;
    if (options::value_of<options::debug_level>(options)>0)
      std::cout << "Degree of rational AAA polynomial " << PEP.num_terms() << std::endl ;

    eigs::info information ;
    auto linearization = CORK::linearization::make_cork_linearization( PEP, information ) ;

    // Describe which eigenvalues are wanted.
    //
    typedef decltype(std::abs(typename EigDomain::value_type()))                                              real_type ;
    CORK::eigenvalue_selector::inside_domain< typename std::decay<decltype(eig_domain)>::type > eigenvalue_selector( eig_domain, options::value_of<options::krylov_breakdown_tolerance<real_type>>(options) ) ;

    auto shift_selector = produce_shift_selection(eigenvalue_selector,options);
    auto const& shifts = shift_selector.shifts() ;
    shift_generator::multiple_shifts< decltype(shifts) > shift_generator( shifts, shift_selector.multiplicity() ) ;

    typedef typename NEP::template value_type_for< typename std::decay<decltype(shift_generator)>::type::value_type> value_type ;

    // Create the CORK triple.
    CORK::krylov::cork_quadruple<real_type,decltype(options)> quadruple( linearization, n_wanted, options ) ;

    // Choose a random initial vector of full rank of the form
    // Q U^T where Q in R^{n\times 2} and U in R^2
    CORK::krylov::random_initial_vector initial_vector( linearization.size_of_basis() ) ;

    // Compute the Schur decomposition. Stop when at least n_wanted eigenvalues have converged.
    eigs::stop_criterion<value_type,NEP,options::stop_criterion<real_type>> stop_criterion( problem, norm_est, CORK::options::value_of<options::stop_criterion<real_type>>(options) ) ;
    CORK::eigs::cork( linearization, eigenvalue_selector, shift_generator, initial_vector, quadruple, stop_criterion, options ) ;

    // Compute the eigenvalues and eigenvectors

    auto [eigenvalues, eigenvectors, resid ] = CORK::eigs::cork_eigen_pairs( problem, std::tuple(quadruple,linearization,eigenvalue_selector,&information) ) ;

    return std::tuple(eigenvalues, eigenvectors, resid,information) ;
  } // eigenvalues_and_vectors()

  template <typename NEP, typename Domain, typename Options>
  decltype (auto) compute_eigenvalues_and_vectors( NEP const& problem, Domain const& aaa_domain, int n_wanted, Options const& options ) {
    return compute_eigenvalues_and_vectors( problem, aaa_domain, aaa_domain, n_wanted, options ) ;
  }

} } } // namespace CORK::eigs::AAA_CORK

#endif
