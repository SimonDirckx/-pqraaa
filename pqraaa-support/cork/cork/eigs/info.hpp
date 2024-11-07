//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_eigs_info_hpp
#define cork_eigs_info_hpp

#include <cork/krylov/info.hpp>
#include <cork/approximation/info.hpp>
#include <cork/krylov/timings.hpp>
#include <limits>

namespace CORK { namespace eigs {

  struct info
  : krylov::info
  {
    inline info()
    // number of converged eigenvalues
    : number_converged(0)
    // number of converged eigenvalues that the user is interested in ????
    , number_converged_and_wanted(0)
    , number_of_restarts(0)
    , total_time(0.)
    , minimum_subspace_dimension(-1)
    , maximum_subspace_dimension(-1)
    , recurrence_error( std::numeric_limits<double>::quiet_NaN() )
    , orthogonalization_error( recurrence_error )
    {}

    int                         number_converged ;
    int                         number_converged_and_wanted ;
    int                         number_of_restarts ;
    double                      total_time ;
    int                         minimum_subspace_dimension ;
    int                         maximum_subspace_dimension ;
    double                      recurrence_error ;
    double                      orthogonalization_error ;
    approximation::info<double> approximation ;
  } ;

  std::ostream& operator<<( std::ostream& os, info const& inf ) {
    os << static_cast<krylov::info const&>(inf) << "\n" ;
    os << "  number of converged eigenvalues = " << inf.number_converged << "\n" ;
    os << "  number of converged and wanted eigenvalues = " << inf.number_converged_and_wanted << "\n" ;
    os << "  number of restarts = " << inf.number_of_restarts << "\n" ;
    os << "  total time = " << inf.total_time << "\n" ;
    os << "  minimum subspace dimension = " << inf.minimum_subspace_dimension << "\n" ;
    os << "  maximum subspace dimension = " << inf.maximum_subspace_dimension << "\n" ;
    os << "  recurrence error = " << inf.recurrence_error << "\n" ;
    os << "  orthogonalization error = " << inf.orthogonalization_error << "\n" ;
    os << inf.approximation ;
    return os ;
  }
   
} } // namespace CORK::eigs


#endif
