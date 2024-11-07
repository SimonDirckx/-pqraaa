//  (C) Copyright Roel Van Beeumen & Karl Meerbergen 2019.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_exception_option_out_of_bounds_hpp
#define cork_exception_option_out_of_bounds_hpp

#include <exception>
#include <string>

namespace CORK { namespace exception {

struct option_out_of_bounds
: std::runtime_error
{
  option_out_of_bounds( std::string const& option, std::string const& info )
  : std::runtime_error( std::string("CORK: Options ") + option + std::string(" is out of bounds: ") + info )
  {}
} ;

} } // namespace CORK::exception

#endif
