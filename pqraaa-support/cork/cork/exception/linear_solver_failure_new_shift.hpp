//  (C) Copyright & Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_exception_linear_solver_failure_new_shift_hpp
#define cork_exception_linear_solver_failure_new_shift_hpp

#include <exception>

namespace CORK { namespace exception {

template <typename T>
class linear_solver_failure_new_shift
: std::runtime_error
{
  public:
    linear_solver_failure_new_shift( T const& shift )
    : std::runtime_error( "CORK: Failure of the linear system solver" )
    , shift_( shift )
    {}

  public:
    T const& shift() const { return shift_ ; }

  private:
    T shift_ ;
} ;

} } // namespace CORK::exception

#endif
