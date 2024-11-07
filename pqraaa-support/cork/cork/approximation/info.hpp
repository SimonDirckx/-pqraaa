//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_approximation_info_hpp
#define cork_approximation_info_hpp

#include <limits>
#include <iostream>
#include <vector>

namespace CORK { namespace approximation {

  template <typename T>
  struct info
  {
    inline info(std::string const& name="NO INFORMATION")
    : approximation(name)
    , size(0)
    , error(std::numeric_limits<float>::infinity() )
    {}

    typedef T                  value_type ;
    std::string                approximation ;
    int                        size ;
    value_type                 error ;
    std::vector< std::string > warnings ;
  } ;

  template <typename T>
  inline std::ostream& operator<<( std::ostream& os, info<T> const& inf ) {
    os << "  number of terms of " << inf.approximation << " : " << inf.size << "\n" ;
    os << "  error of " << inf.approximation << " : " << inf.error << "\n" ;
    os << "  warnings for " << inf.approximation << " : " ;
    for (auto it=inf.warnings.begin(); it!=inf.warnings.end(); ++it) os << "{" << *it << "}" ;
    return os ;
  }
   
} } // namespace CORK::approximation


#endif
