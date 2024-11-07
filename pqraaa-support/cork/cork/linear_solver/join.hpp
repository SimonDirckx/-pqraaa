//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_linear_solver_linear_solver_join_hpp
#define cork_linear_solver_linear_solver_join_hpp

#include <cork/linear_solver/linear_solver.hpp>
#include <cork/coefficient_matrices/join.hpp>
#include <cork/linear_solver/dense.hpp>
#include <cork/linear_solver/sparse.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace linear_solver {

  template <typename T, typename Front>
  struct linear_solver_traits< T, coefficient_matrices::join<Front, Front> >
  : CORK::linear_solver::linear_solver_traits<T, Front>
  {
     static auto apply( coefficient_matrices::join<Front,Front> const& join ) {
       return CORK::linear_solver::linear_solver_traits<T, Front>::apply( join.front() ) ;
     }
  } ;

  template <typename T, typename Front, typename Tail>
  struct linear_solver_traits< T, coefficient_matrices::join<Front, Tail>
                             , typename std::enable_if< glas2::is<glas2::CoordinateSparseMatrix, typename std::decay<Front>::type::matrix_type>::value
                                                     && glas2::is<glas2::DenseMatrix, typename std::decay<Tail>::type::matrix_type>::value
                                                     >::type >
  {
    typedef typename std::decay<Front>::type front_type ;
    typedef typename std::decay<Tail>::type  tail_type ;
    typedef typename std::common_type< T, typename std::common_type< typename front_type::value_type
                                                                   , typename tail_type::value_type
                                                                   >::type >::type       value_type ;
    typedef CORK::linear_solver::dense< value_type >                                     type ;

    template <typename J>
    static type apply( J const& join ) {
      return type( join.num_rows() ) ;
    }
  } ;

  template <typename T, typename Front, typename Tail>
  struct linear_solver_traits< T, coefficient_matrices::join<Front, Tail>
                             , typename std::enable_if< glas2::is<glas2::CoordinateSparseMatrix, typename std::decay<Tail>::type::matrix_type>::value
                                                     && glas2::is<glas2::DenseMatrix, typename std::decay<Front>::type::matrix_type>::value
                                                     >::type >
  : linear_solver_traits< T, coefficient_matrices::join<Tail, Front> >
  {} ;

} } // namespace CORK::linear_solver

#endif
