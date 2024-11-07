//  (C) Copyright Roel Van Beeumen & Karl Meerbergen 2019.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_exception_initial_vector_hpp
#define cork_exception_initial_vector_hpp

#include <exception>
#include <string>

namespace CORK { namespace exception {

struct initial_vector
: std::runtime_error
{
  inline initial_vector()
  : std::runtime_error( std::string("CORK: Initial vector is zero") )
  {}
} ;

} } // namespace CORK::exception

#endif
