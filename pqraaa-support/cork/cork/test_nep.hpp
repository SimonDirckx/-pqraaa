//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_test_nep_hpp
#define cork_test_nep_hpp

#include <cork/linear_solver/user_defined.hpp>
#include <cork/matrix_valued_function/linear_solver.hpp>
#include <cork/linear_solver/user_defined.hpp>
#include <cork/utility/ref.hpp>
#include <cork/utility/value_type.hpp>
#include <cork/concept/stl_sequence_container.hpp>
#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>

namespace CORK { 

 /**
  * Test the consistency of the multiply_add and solve functions of the Nonlinear Eigenvalue Problem.
  * 
  * @param problem        The nonlinear eigenvalue problem, which contains a multiply_add and a solve function
  * @param lambda_begin   The beginning of a forward iterator
  * @param lambda_end     The end of this forward iterator
  */
  template <typename NEP, typename TP>
#ifdef CORK_USE_CONCEPTS
  requires CORK::STLForwardIterator<TP>
#endif
  void test_nep( NEP problem, TP lambda_begin, TP lambda_end ) {
    typedef typename CORK::deref_type<NEP>::type::template value_type_for< typename std::decay<decltype(*lambda_begin)>::type > value_type ;
    glas2::vector< value_type > x( CORK::deref(problem).size() ) ;
    glas2::vector< value_type > y( CORK::deref(problem).size() ) ;

    std::cout << typeid(std::ref(CORK::deref(problem))).name() << std::endl;

    auto solver = CORK::matrix_valued_function::make_linear_solver<value_type>( &CORK::deref(problem) ) ;

    for (auto lambda_p=lambda_begin; lambda_p!=lambda_end; ++lambda_p) {
      glas2::randomize(x) ;
      fill(x,1.);
      y = -x ;
      solver.prepare_solve( *lambda_p ) ;
      solver.solve( y ) ;
      decltype(std::abs(value_type())) norm_x = norm_2( x ) ;
      CORK::deref(problem).multiply_add( *lambda_p, y, x ) ;
      std::cout << "Norm of error for lambda = " << *lambda_p << " is " << norm_2(x) << " / " << norm_x << std::endl ;
    }
  } // test_nep()

} // namespace CORK

#endif
