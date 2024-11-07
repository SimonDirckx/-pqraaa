//  (C) Copyright Roel Van Beeumen & Karl Meerbergen 2019.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_exception_not_implemented_hpp
#define cork_exception_not_implemented_hpp

#include <exception>

namespace CORK { namespace exception {

struct not_implemented
: std::runtime_error
{
  not_implemented( std::string const& str = "This option is not implemented yet." )
  : std::runtime_error( "CORK: " + str )
  {}
} ;

} } // namespace CORK::exception

#endif
