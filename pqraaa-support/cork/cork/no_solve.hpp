//  (C) Copyright Karl Meerbergen & Dries De Samblanx, 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_no_solve
#define cork_no_solve

#include <cork/vector.hpp>
#include <cork/exception/not_implemented.hpp>

namespace CORK {

  struct no_solve {
    template <typename T>
    void operator() ( T const& shift, CORK::vector<T> , bool ) const {
      throw exception::not_implemented( "a branch of CORK is called requiring a linear solver" ) ;
    }
  } ;

} // namespace CORK

#endif

