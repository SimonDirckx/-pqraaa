//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_matrix_iterator_newton_hpp
#define cork_matrix_iterator_newton_hpp

#include <cork/matrix_iterator/matrix_iterator.hpp>
#include <cork/basis/newton.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace matrix_iterator {

  // MatrixHandler contains the matrix coefficients:
  //
  // A_i i=0,...,grade-1 contains the A-matrices
  // A_d contains B_{d-1}
  template <typename Points, typename I>
  class matrix_iterator< basis::newton< Points, I > >
  {
    public:
      typedef basis::newton< Points, I >      basis_type ;
      typedef typename basis_type::size_type  grade_type ;

    private:
      typedef typename basis_type::points_type points_type ;

    public:
      matrix_iterator( basis::newton< Points, I > const& basis )
      : basis_( basis )
      {}

    public:
      grade_type num_CORK_coefs() const {
        return basis_.num_terms()-1 ;
      }

      grade_type num_coefs() const {
        return basis_.num_terms() ;
      }

    public:
      template <typename CoefficientMatrices, typename ZLambda>
      void schedule_a_0( CoefficientMatrices& accumulator, ZLambda const& z_lambda ) const {
        for (int i=0; i<num_CORK_coefs(); ++i) {
          accumulator.schedule( i, z_lambda(i) ) ;
        }
        accumulator.schedule( num_CORK_coefs(), -basis_.points()(num_CORK_coefs()-1) * z_lambda(num_CORK_coefs()-1) ) ;
      }

      template <typename CoefficientMatrices, typename ZLambda>
      void schedule_a_1( CoefficientMatrices& accumulator, ZLambda const& z_lambda ) const {
        for (int i=1; i<num_CORK_coefs(); ++i) {
          accumulator.schedule( i, z_lambda(i) ) ;
        }
        accumulator.schedule( num_CORK_coefs(), -basis_.points()(num_CORK_coefs()-1) * z_lambda(num_CORK_coefs()-1) ) ;
      }

      template <typename CoefficientMatrices, typename ZLambda>
      void schedule_b_0( CoefficientMatrices& accumulator, ZLambda const& z_lambda ) const {
        accumulator.schedule( num_CORK_coefs(), -z_lambda(num_CORK_coefs()-1) ) ;
      }

      template <typename CoefficientMatrices, typename ZLambda>
      void schedule_b_1( CoefficientMatrices& accumulator, ZLambda const& z_lambda ) const {
        accumulator.schedule( num_CORK_coefs(), -z_lambda(num_CORK_coefs()-1) ) ;
      }

    public:
      basis_type        basis_ ;
  } ; // class newton


} } // namespace CORK::matrix_iterator

#endif
