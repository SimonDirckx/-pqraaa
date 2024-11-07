//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_linear_solver_linear_solver_combined_hpp
#define cork_linear_solver_linear_solver_combined_hpp

#include <cork/linear_solver/linear_solver.hpp>
#include <cork/coefficient_matrices/combined.hpp>
#include <cork/linear_solver/combined.hpp>
#include <type_traits>

namespace CORK { namespace linear_solver {

  template <typename T, typename CoefficientMatrices, typename Combinations>
  struct linear_solver_traits< T, coefficient_matrices::combined<CoefficientMatrices, Combinations> >
  {
    typedef typename std::common_type< T, typename std::decay<Combinations>::type::value_type>::type value_type ;

    typedef typename linear_solver_traits< value_type, typename std::decay<CoefficientMatrices>::type >::type type ;

    static auto apply( coefficient_matrices::combined<CoefficientMatrices,Combinations> const& combined ) {
      return make_linear_solver<value_type>( combined.coefficient_matrices() ) ;
    }
  } ;

} } // namespace CORK::linear_solver

#endif
