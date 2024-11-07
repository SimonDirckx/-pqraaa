//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompany_glasing file 
//  LICENSE_1_0.txt)

#ifndef cork_linear_solver_any_glas_hpp
#define cork_linear_solver_any_glas_hpp

#include <cork/coefficient_matrices/any_glas.hpp>
#include <cork/linear_solver/linear_solver.hpp>
#include <cork/linear_solver/dense.hpp>
#include <cork/linear_solver/sparse.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace linear_solver {

  template <typename T, typename Sequence>
  struct linear_solver_traits< T, coefficient_matrices::any_glas<Sequence>
                             , typename std::enable_if< glas2::is<glas2::DenseMatrix, typename coefficient_matrices::any_glas<Sequence>::matrix_type>::value >::type
                             >
  {
    typedef typename coefficient_matrices::any_glas<Sequence>::matrix_type                      matrix_type ;
    typedef typename std::common_type< T, typename matrix_type::value_type>::type               value_type ;
    typedef CORK::linear_solver::dense< value_type >                                            type ;

    template <typename J>
    static type apply( J const& any_glas ) {
      return type( CORK::deref(any_glas).num_rows() ) ;
    }
  } ;

  template <typename T, typename Sequence>
  struct linear_solver_traits< T, coefficient_matrices::any_glas<Sequence>
                             , typename std::enable_if< glas2::is<glas2::CoordinateSparseMatrix, typename coefficient_matrices::any_glas<Sequence>::matrix_type>::value >::type
                             >
  {
    typedef typename coefficient_matrices::any_glas<Sequence>::matrix_type                      matrix_type ;
    typedef typename std::common_type< T, typename matrix_type::value_type>::type               value_type ;
    typedef CORK::linear_solver::sparse< value_type >                                           type ;

    template <typename J>
    static type apply( J const& any_glas ) {
      return type( CORK::deref(any_glas).num_rows() ) ;
    }
  } ;

  template <typename T, typename Sequence>
  struct linear_solver_traits< T, coefficient_matrices::any_glas<Sequence>
                             , typename std::enable_if< glas2::is<glas2::CompressedSparseMatrix, typename coefficient_matrices::any_glas<Sequence>::matrix_type>::value >::type
                             >
  {
    typedef typename coefficient_matrices::any_glas<Sequence>::matrix_type                      matrix_type ;
    typedef typename std::common_type< T, typename matrix_type::value_type>::type               value_type ;
    typedef CORK::linear_solver::sparse< value_type >                                           type ;

    template <typename J>
    static type apply( J const& any_glas ) {
      return type( CORK::deref(any_glas).num_rows() ) ;
    }
  } ;

} } // namespace CORK::linear_solver

#endif
