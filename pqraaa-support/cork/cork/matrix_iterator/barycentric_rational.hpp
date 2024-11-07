//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_matrix_iterator_barycentric_rational_hpp
#define cork_matrix_iterator_barycentric_rational_hpp

#include <cork/matrix_iterator/matrix_iterator.hpp>
#include <cork/basis/barycentric_rational.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace matrix_iterator {

  // MatrixHandler contains the matrix coefficients:
  //
  // A_i i=0,...,grade-1 contains the A-matrices
  // A_d contains B_{d-1}
  template <typename Weights, typename Poles, typename I>
  class matrix_iterator< basis::barycentric_rational< Weights, Poles, I > >
  {
    public:
      typedef basis::barycentric_rational< Weights, Poles, I >      basis_type ;
      typedef typename basis_type::size_type  grade_type ;

    private:
      typedef typename basis_type::weights_type weights_type ;

    public:
      matrix_iterator( basis::barycentric_rational< Weights, Poles, I > const& basis )
      : basis_( basis )
      {}

    public:
      grade_type num_CORK_coefs() const {
        return basis_.num_terms() ;
      }

      grade_type num_coefs() const {
        return basis_.num_terms() ;
      }

    public:
      template <typename CoefficientMatrices4CORK, typename ZLambda>
      void schedule_a_0( CoefficientMatrices4CORK& accumulator, ZLambda const& z_lambda ) const {
        // 1st block is zero, so start from i=1
        accumulator.schedule( 0, z_lambda(0) ) ;
        for (int i=1; i<num_CORK_coefs(); ++i) {
#ifndef CORK_BARYCENTRIC_RATIONAL_WEIGHTED
          accumulator.schedule( i, z_lambda(i) * basis_.weights()(i-1) ) ;
#else
          accumulator.schedule( i, z_lambda(i) ) ;
#endif
        }
      }

      template <typename CoefficientMatrices4CORK, typename ZLambda>
      void schedule_a_1( CoefficientMatrices4CORK& accumulator, ZLambda const& z_lambda ) const {
        for (int i=1; i<num_CORK_coefs(); ++i) {
#ifndef CORK_BARYCENTRIC_RATIONAL_WEIGHTED
          accumulator.schedule( i, z_lambda(i) * basis_.weights()(i-1) ) ;
#else
          accumulator.schedule( i, z_lambda(i) ) ;
#endif
        }
      }

      template <typename CoefficientMatrices4CORK, typename ZLambda>
      void schedule_b_0( CoefficientMatrices4CORK& accumulator, ZLambda const& z_lambda ) const { }

      template <typename CoefficientMatrices4CORK, typename ZLambda>
      void schedule_b_1( CoefficientMatrices4CORK& accumulator, ZLambda const& z_lambda ) const { }

    public:
      basis_type        basis_ ;
  } ; // class barycentric_rational


} } // namespace CORK::matrix_iterator

#endif
