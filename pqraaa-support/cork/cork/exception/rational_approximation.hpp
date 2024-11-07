//  (C) Copyright Roel Van Beeumen & Karl Meerbergen 2019.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_exception_rational_approximation_hpp
#define cork_exception_rational_approximation_hpp

#include <exception>
#include <string>

namespace CORK { namespace exception {

struct rational_approximation
: std::runtime_error
{
  rational_approximation( std::string const& s )
  : std::runtime_error( "CORK: " + s )
  {}
} ;

} } // namespace CORK::exception

#endif
