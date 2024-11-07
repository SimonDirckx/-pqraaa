//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_eigs_pep_cork_hpp
#define cork_eigs_pep_cork_hpp

#include <cork/eigs/cork.hpp>
#include <cork/eigs/options.hpp>
#include <cork/eigs/info.hpp>
#include <cork/krylov/random_initial_vector.hpp>
#include <cork/shift_generator/multiple_poles.hpp>
#include <cork/linearization/cork_linearization.hpp>
#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <iosfwd>

namespace CORK { namespace eigs { namespace PEP_CORK {

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


  template <typename T=double>
  struct options {
    typedef decltype(std::abs(T())) real_type ;

    options()
    : number_of_shifts(2)
    , number_of_shift_changes_per_restart(1)
    , given_shifts(0)
    {}

    eigs::options<real_type>                     krylov_options ;
    int                                          number_of_shifts ;
    int                                          number_of_shift_changes_per_restart ;
    glas2::vector<T>                             given_shifts ;
  } ;


  template <typename PEP, typename Geometry, typename Options, typename RealWanted=std::false_type >
  decltype (auto) compute_eigenvalues_and_vectors( PEP const& problem, Geometry const& eigenvalue_selector, int n_wanted, Options const& options, RealWanted=std::false_type() ) {
    typedef typename PEP::template value_type< typename Geometry::value_type > value_type ;

/*    typedef typename std::conditional< RealWanted::value && !std::is_same< typename PEP::template value_type< decltype(std::abs(value_type())) >, value_type >::value
                                     , typename PEP::template value_type< decltype(std::abs(value_type())) >
                                     , value_type
                                     >::type q_value_type ;
                                     */
    typedef value_type q_value_type ;

    std::vector< std::string > warnings ;

    int number_of_shifts = options.number_of_shifts ;
    if (number_of_shifts<1) {
      number_of_shifts = 1 ;
      std::cout <<"CORK WARNING: Options variable number_of_shifts is too small. Used value 1." << std::endl ;
      warnings.push_back( "Options variable number_of_shifts is smaller than 1." ) ;
    }

    int number_of_shift_changes_per_restart = options.number_of_shift_changes_per_restart ;
    if (number_of_shift_changes_per_restart<1) {
      number_of_shift_changes_per_restart = 1 ;
      std::cout <<"CORK WARNING: Options variable number_of_shift_changes_per_restart is too small. Used value 1." << std::endl ;
      warnings.push_back( "Options variable number_of_shift_changes_per_restart is smaller than 1." ) ;
    }

    eigs::info information ;
    auto linearization = CORK::linearization::make_cork_linearization( problem, information ) ;

    if (options.krylov_options.debug_level>0)
      std::cout << "shifts to be used " << eigenvalue_selector.poles() << std::endl ;

    // Create the CORK triple.
    CORK::krylov::cork_quadruple<q_value_type,decltype(options.krylov_options)> quadruple( linearization, eigenvalue_selector, options.krylov_options ) ;

    // Choose a random initial vector of full rank of the form
    // Q U^T where Q in R^{n\times 2} and U in R^2
    CORK::krylov::random_initial_vector initial_vector( linearization.size_of_basis() ) ;

    // Choose poles in blocks
    CORK::shift_generator::multiple_poles< decltype(eigenvalue_selector.poles()) const& > shift_generator( eigenvalue_selector.poles(), std::max<int>(1,quadruple.k_max() / (number_of_shift_changes_per_restart/**eigenvalue_selector.poles().size()*/)) ) ;

    // Compute the Schur decomposition.
    CORK::eigs::cork( linearization, eigenvalue_selector, shift_generator, initial_vector, quadruple, options.krylov_options ) ;

    // Compute the eigenvalues and eigenvectors
    eigen_system< value_type > eigen( problem.size(), information.number_converged_and_wanted ) ;
    eigen.shifts.resize( eigenvalue_selector.poles().size() ) ;
    eigen.shifts = eigenvalue_selector.poles() ;

    CORK::eigs::cork_eigen( linearization, eigenvalue_selector, quadruple, eigen.eigen_values, eigen.eigen_vectors ) ;

    eigen.information = information ;
    return eigen ;
  } // eigenvalues_and_vectors()

} } } // namespace CORK::eigs::PEP_CORK

#endif
