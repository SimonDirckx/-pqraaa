//  (C) Copyright Roel Van Beeumen & Karl Meerbergen 2019.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_exception_linear_solver__failure_hpp
#define cork_exception_linear_solver__failure_hpp

#include <stdexcept>

namespace CORK { namespace exception {

struct linear_solver_failure
: std::runtime_error
{
  linear_solver_failure()
  : std::runtime_error( "CORK: Failure of the linear system solver" )
  {}
} ;

} } // namespace CORK::exception

#endif
