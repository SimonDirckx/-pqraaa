//  (C) Copyright Karl Meerbergen 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_matrix_iterator_accumulated_partial_fractions_real_hpp
#define cork_matrix_iterator_accumulated_partial_fractions_real_hpp

#include <cork/matrix_iterator/matrix_iterator.hpp>
#include <cork/basis/accumulated_partial_fractions_real.hpp>
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
  class matrix_iterator< basis::accumulated_partial_fractions_real< Weights, Poles, I > >
  {
    public:
      typedef basis::accumulated_partial_fractions_real< Weights, Poles, I >      basis_type ;
      typedef typename basis_type::size_type  grade_type ;

    public:
      template <typename ValueType>
      using value_type = ValueType ;

    public:
      matrix_iterator( basis_type const& basis )
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
      template <typename CoefficientMatrices, typename ZLambda>
      void schedule_a_0( CoefficientMatrices& iterator, ZLambda const& z_lambda ) const {
        for (int i=0; i<num_CORK_coefs(); ++i) {
          iterator.schedule( i, z_lambda(i) ) ;
        }
      }

      template <typename CoefficientMatrices, typename ZLambda>
      void schedule_a_1( CoefficientMatrices& iterator, ZLambda const& z_lambda ) const {
        for (int i=1; i<num_CORK_coefs(); ++i) {
          iterator.schedule( i, z_lambda(i) ) ;
        }
      }

      template <typename CoefficientMatrices, typename ZLambda>
      void schedule_b_0( CoefficientMatrices&, ZLambda const& z_lambda ) const { }

      template <typename CoefficientMatrices, typename ZLambda>
      void schedule_b_1( CoefficientMatrices&, ZLambda const& z_lambda ) const { }

    public:
      basis_type        basis_ ;
  } ; // class partial_fractions


} } // namespace CORK::matrix_iterator

#endif
