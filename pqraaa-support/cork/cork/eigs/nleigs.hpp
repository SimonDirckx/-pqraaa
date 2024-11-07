//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_eigs_nleigs_hpp
#define cork_eigs_nleigs_hpp

#include <cork/eigs/cork.hpp>
#include <cork/eigs/options.hpp>
#include <cork/eigs/info.hpp>
#include <cork/krylov/random_initial_vector.hpp>
#include <cork/linearization/union_of_bases.hpp>
#include <cork/linearization/cork_linearization.hpp>
#include <cork/linearization/rational_newton.hpp>
#include <cork/linearization/monomial.hpp>
#include <cork/linearization/combined.hpp>
#include <cork/matrix_valued_function/potential_theory.hpp>
#include <cork/eigenvalue_selector/inside_domain.hpp>
#include <cork/shift_generator/greedy.hpp>
#include <cork/domain/centers.hpp>
#include <glas2/bindings/vector.hpp>
#include <boost/numeric/bindings/std/vector.hpp>

namespace CORK { namespace eigs { namespace NLEIGS {

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


  // Options should contain which NLEIGS, which shift_generator to choose
  template <typename Scenario, typename T=typename Scenario::value_type>
  struct options {
    typedef decltype(std::abs(T())) real_type ;

    options( Scenario const& shift_generator )
    : shift_generator( shift_generator )
    {}

    matrix_valued_function::potential_theory_options<real_type> potential_theory_options ;
    eigs::options<real_type>                       krylov_options ;
    Scenario                                       shift_generator ;
  } ;

  template <typename Scenario>
  decltype (auto) make_options( Scenario const& shift_generator ) { return options<Scenario>( shift_generator ) ; }


  template <typename NEP, typename Sigma, typename Xi, typename Options>
  decltype (auto) compute_eigenvalues_and_vectors( NEP const& problem, Sigma const& sigma, Xi const& xi, int n_wanted, Options const& options ) {
    return compute_eigenvalues_and_vectors( problem, sigma, xi, sigma, n_wanted, options ) ;
  }

  template <typename NEP, typename Sigma, typename Xi, typename EigDomain, typename Options>
  decltype (auto) compute_eigenvalues_and_vectors( NEP const& problem, Sigma const& sigma, Xi const& xi, EigDomain const& eig_domain, int n_wanted, Options const& options ) {
    typedef typename NEP::template value_type< typename EigDomain::value_type > value_type ;

    std::vector< std::string > warnings ;

    auto PEP = CORK::matrix_valued_function::potential_theory< value_type >( problem, sigma, xi, options.potential_theory_options ) ;
    if (options.krylov_options.debug_level>0)
      std::cout << "Degree of rational NLEIGS polynomial " << PEP.num_terms() << std::endl ;

    eigs::info information ;
    auto linearization = CORK::linearization::make_cork_linearization( PEP, information ) ;

    // Describe which eigenvalues are wanted.
    //
    CORK::eigenvalue_selector::inside_domain< typename std::decay<decltype(eig_domain)>::type > eigenvalue_selector( eig_domain, n_wanted, options.krylov_options.relative_tolerance ) ;

    // Create the CORK triple.
    CORK::krylov::cork_quadruple<value_type,decltype(options.krylov_options)> quadruple( linearization, eigenvalue_selector, options.krylov_options ) ;

    // Choose a random initial vector of full rank of the form
    // Q U^T where Q in R^{n\times 2} and U in R^2
    CORK::krylov::random_initial_vector initial_vector( linearization.size_of_basis() ) ;

    // Compute the Schur decomposition. Stop when at least 5 eigenvalues have converged.
    CORK::eigs::cork( linearization, eigenvalue_selector, options.shift_generator, initial_vector, quadruple, options.krylov_options ) ;

    // Compute the eigenvalues and eigenvectors
    eigen_system< value_type > eigen( PEP.size(), information.number_converged_and_wanted ) ;
    eigen.shifts.resize( options.shift_generator.used_shifts().size() ) ;
    eigen.shifts = glas2::bind_vector( options.shift_generator.used_shifts() ) ;

    CORK::eigs::cork_eigen( linearization, eigenvalue_selector, quadruple, eigen.eigen_values, eigen.eigen_vectors ) ;

    eigen.information = information ;
    return eigen ;
  } // eigenvalues_and_vectors()

} } } // namespace CORK::eigs::NLEIGS

#endif
