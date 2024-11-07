//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_USE_PEP_LINEAR_SOLVER_HPP
#define CORK_OPTIONS_USE_PEP_LINEAR_SOLVER_HPP

#include<limits>
#include<tuple>

namespace CORK { namespace options {

  class use_pep_linear_solver {
    public:
      use_pep_linear_solver()
      {}

      template <typename ...Ts>
      use_pep_linear_solver( std::tuple<Ts...> const& options)
      {}

      use_pep_linear_solver( bool v )
      {}

      bool value() const { return std::true_type() ; }
  } ;

} } // CORK::options

#endif
