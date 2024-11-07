//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_linearization_info_hpp
#define cork_linearization_info_hpp

#include <limits>
#include <cmath>
#include <string>
#include <iostream>
#include <vector>

namespace CORK { namespace linearization {

  struct info {
    inline info()
    : number_of_solves(0)
    , number_of_factorizations(0)
    , time_of_solves(0.)
    , time_of_factorizations(0.)
    , time_of_matvecs(0.)
    , linearization("NO INFORMATION")
    , error("")
    {}

    int                        number_of_solves ;
    int                        number_of_factorizations ;
    double                     time_of_solves ;
    double                     time_of_factorizations ;
    double                     time_of_matvecs ;
    std::string                linearization ;
    std::string                error ;
    std::vector< std::string > warnings ;
  } ;

  inline std::ostream& operator<<( std::ostream& os, info const& inf ) {
    os << "CORK information:\n";
    os << "  type of linearization = " << inf.linearization << "\n" ;
    os << "  number of linear solves = " << inf.number_of_solves << "\n" ;
    os << "  number of factorizations = " << inf.number_of_factorizations << "\n" ;
    os << "  time for linear solves = " << inf.time_of_solves << "\n" ;
    os << "  time for factorizations = " << inf.time_of_factorizations << "\n" ;
    os << "  time for matrix vector products = " << inf.time_of_matvecs << "\n" ;
    os << "  error message: " << inf.error << "\n" ;
    os << "  warnings: " ;
    for (auto it=inf.warnings.begin(); it!=inf.warnings.end(); ++it) os << *it << ", " ;
    return os ;
  }
   
} } // namespace CORK::linearization


#endif
