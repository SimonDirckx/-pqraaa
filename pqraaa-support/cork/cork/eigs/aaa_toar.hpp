//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_eigs_aaa_toar_hpp
#define cork_eigs_aaa_toar_hpp

#include <cork/eigs/info.hpp>
#include <cork/eigs/options.hpp>
#include <cork/eigs/toar.hpp>
#include <cork/krylov/random_initial_vector.hpp>
#include <cork/linearization/union_of_bases.hpp>
#include <cork/linearization/cork_linearization.hpp>
#include <cork/linearization/barycentric_rational.hpp>
#include <cork/linearization/monomial.hpp>
#include <cork/linearization/combined.hpp>
#include <cork/matrix_valued_function/sv_aaa.hpp>
#include <cork/eigenvalue_selector/inside_domain_given_shifts.hpp>
#include <cork/domain/center.hpp>
#include <cork/approximation/sv_aaa.hpp>

namespace CORK { namespace eigs { namespace AAA_TOAR {

  template <typename T>
  struct eigen_system {
    eigen_system( int n, int n_cvgd )
    : eigen_values(n_cvgd)
    , eigen_vectors(n,n_cvgd)
    {}

    glas2::vector<T> eigen_values ;
    glas2::matrix<T> eigen_vectors ;
    eigs::info       information ;
    T                shift ;
  } ;


  template <typename T=double>
  struct options {
    typedef decltype(std::abs(T())) real_type ;
    matrix_valued_function::aaa_options<real_type > aaa_options ;
    eigs::options<real_type>                        krylov_options ;
  } ;


  template <typename NEP, typename Domain, typename Options>
  decltype (auto) compute_eigenvalues_and_vectors( NEP const& problem, Domain const& domain, int n_wanted, Options const& options ) {
    typedef typename std::common_type< typename NEP::value_type, typename Domain::value_type>::type value_type ;

    auto PEP = CORK::matrix_valued_function::sv_aaa< value_type >( problem, domain, options.aaa_options ) ;
    if (options.krylov_options.debug_level>0)
      std::cout << "Degree of rational AAA polynomial " << PEP.num_terms() << std::endl ;

    CORK::eigs::info information ;
    auto linearization = CORK::linearization::make_cork_linearization< value_type >( PEP, information ) ;

    // Describe which eigenvalues are wanted.
    glas2::vector< value_type > shifts(1) ; shifts(0) = domain::center( domain ) ;
    CORK::eigenvalue_selector::inside_domain_given_shifts< decltype(shifts), typename std::decay<decltype(domain)>::type > eigenvalue_selector( shifts, domain, n_wanted, options.krylov_options.relative_tolerance ) ;

    if (options.krylov_options.debug_level>0) std::cout << "Krylov shift is " << shifts(0) << std::endl ;

    // Create the CORK triple.
    CORK::krylov::toar_triple<value_type> triple( linearization, eigenvalue_selector, options.krylov_options ) ;

    // Choose a random initial vector of full rank of the form
    // Q U^T where Q in R^{n\times 2} and U in R^2
    CORK::krylov::random_initial_vector initial_vector( linearization.size_of_basis() ) ;

    // Compute the Schur decomposition. Stop when at least 5 eigenvalues have converged.
    CORK::eigs::toar( linearization, eigenvalue_selector, initial_vector, triple, options.krylov_options ) ;

    // Compute the eigenvalues and eigenvectors
    eigen_system< value_type > eigen( PEP.size(), information.number_converged_and_wanted ) ;
    eigen.information = information ;
    eigen.shift = shifts(0) ;
    CORK::eigs::toar_eigen( linearization, eigenvalue_selector, triple, eigen.eigen_values, eigen.eigen_vectors ) ;

    eigen.information = information ;
    return eigen ;
  } // AAA_TOAR_eigenvalues_and_vectors()

} } } // namespace CORK::eigs::AAA_TOAR

#endif
