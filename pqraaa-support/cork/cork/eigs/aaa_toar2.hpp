//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_eigs_aaa_toar2_hpp
#define cork_eigs_aaa_toar2_hpp

#ifdef cork_eigs_aaa_toar_hpp
#error aaa_toar2_hpp and aaa_toar.hpp cannot be included together
#endif

#include <cork/eigs/toar.hpp>
#include <cork/eigs/options.hpp>
#include <cork/eigs/info.hpp>
#include <cork/krylov/random_initial_vector.hpp>
#include <cork/linearization/union_of_bases.hpp>
#include <cork/linearization/cork_linearization.hpp>
#include <cork/linearization/barycentric_rational.hpp>
#include <cork/linearization/monomial.hpp>
#include <cork/linearization/combined.hpp>
#include <cork/matrix_valued_function/sv_aaa.hpp>
#include <cork/eigenvalue_selector/inside_domain_given_shifts.hpp>
#include <cork/domain/centers.hpp>
#include <cork/approximation/sv_aaa.hpp>
#include <glas2/bindings/vector.hpp>
#include <boost/numeric/bindings/std/vector.hpp>

namespace CORK { namespace eigs { namespace AAA_TOAR {

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


  template <typename T>
  struct options {
    typedef decltype(std::abs(T())) real_type ;

    options()
    {}

    matrix_valued_function::aaa_options<real_type> aaa_options ;
    eigs::options<real_type>                       krylov_options ;
  } ;


  template <typename NEP, typename Domain, typename Options>
  decltype (auto) compute_eigenvalues_and_vectors( NEP const& problem, Domain const& aaa_domain, int n_wanted, Options const& options ) {
    return compute_eigenvalues_and_vectors( problem, aaa_domain, aaa_domain, n_wanted, options ) ;
  }

  template <typename NEP, typename AAADomain, typename EigDomain, typename Options>
  decltype (auto) compute_eigenvalues_and_vectors( NEP const& problem, AAADomain const& aaa_domain, EigDomain const& eig_domain, int n_wanted, Options const& options ) {
    typedef typename NEP::template value_type< typename EigDomain::value_type > value_type ;

    std::vector< std::string > warnings ;

    auto PEP = CORK::matrix_valued_function::sv_aaa< value_type >( problem, aaa_domain, options.aaa_options ) ;
    if (options.krylov_options.debug_level>0)
      std::cout << "Degree of rational AAA polynomial " << PEP.num_terms() << std::endl ;

    eigs::info information ;
    auto linearization = CORK::linearization::make_cork_linearization( PEP, information ) ;

    // Describe which eigenvalues are wanted.
    //
    glas2::vector< value_type > shifts = eig_domain.discretize( 1 ) ;
    
    CORK::eigenvalue_selector::inside_domain_given_shifts< decltype(shifts), typename std::decay<decltype(eig_domain)>::type > eigenvalue_selector( shifts, eig_domain, n_wanted, options.krylov_options.relative_tolerance ) ;

    // Create the CORK triple.
    CORK::krylov::toar_triple<value_type> triple( linearization, eigenvalue_selector, options.krylov_options ) ;

    // Choose a random initial vector of full rank of the form
    // Q U^T where Q in R^{n\times 2} and U in R^2
    CORK::krylov::random_initial_vector initial_vector( linearization.size_of_basis() ) ;

    // Compute the Schur decomposition. Stop when at least 5 eigenvalues have converged.
    CORK::eigs::toar( linearization, eigenvalue_selector, initial_vector, triple, options.krylov_options ) ;

    // Compute the eigenvalues and eigenvectors
    eigen_system< value_type > eigen( PEP.size(), information.number_converged_and_wanted ) ;

    toar_eigen( linearization, eigenvalue_selector, triple, eigen.eigen_values, eigen.eigen_vectors ) ;

    eigen.information = information ;
    return eigen ;
  } // eigenvalues_and_vectors()

} } } // namespace CORK::eigs::AAA_TOAR

#endif
