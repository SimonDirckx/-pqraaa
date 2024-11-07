//  (C) Copyright Roel Van Beeumen & Karl Meerbergen 2019.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_exception_lapack_error_hpp
#define cork_exception_lapack_error_hpp

#include <exception>
#include <cstring>

namespace CORK { namespace exception {

struct lapack_error
: std::runtime_error
{
  lapack_error(std::string const& s)
  : std::runtime_error( std::string("CORK: Failure of a LAPACK routine: ") + s )
  {}
} ;

} } // namespace CORK::exception

#endif
