//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_matrix_polynomial_monomial_matrix_polynomial_hpp
#define cork_matrix_polynomial_monomial_matrix_polynomial_hpp

#include <cork/basis/monomial.hpp>
#include <cork/matrix_iterator/matrix_iterator.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <cmath>
#include <type_traits>

namespace CORK { namespace matrix_iterator {

  // MatrixHandler contains the matrix coefficients:
  //
  // A_i i=0,...,grade-1 contains the A-matrices
  // A_d contains B_{d-1}
  template <typename I>
  class matrix_iterator< basis::monomial<I> >
  {
    public:
      typedef I         grade_type ;

    public:
      template <typename ValueType>
      using value_type = ValueType ;

    public:
      matrix_iterator( basis::monomial<I> const& basis )
      : basis_( basis )
      {
        assert( basis_.num_terms()>0 ) ;
      }

    public:
      grade_type num_coefs() const {
        return basis_.num_terms() ;
      }

    public:
      grade_type num_CORK_coefs() const {
        return std::max(1,basis_.num_terms()-1) ;
      }

    public:
      /*
      template <typename Accumulator, typename ZLambda>
      void schedule_a_0( Accumulator accumulator, ZLambda const& z_lambda ) const {
        assert( accumulator.num_columns()==num_coefs() ) ;
        for (int i=0; i<num_CORK_coefs(); ++i) {
          accumulator( glas2::all(), i ) += z_lambda( i ) ;
        }
      } // schedule_a_0()
      */

      template <typename CoefficientMatrices, typename ZLambda>
      void schedule_a_0( CoefficientMatrices& iterator, ZLambda const& z_lambda ) const {
        for (int i=0; i<num_CORK_coefs(); ++i) {
          iterator.schedule( i, z_lambda( i ) ) ;
        }
      } // schedule_a_0()

      template <typename CoefficientMatrices, typename ZLambda>
      void schedule_a_1( CoefficientMatrices& iterator, ZLambda const& z_lambda ) const {
        for (int i=1; i<num_CORK_coefs(); ++i) {
          iterator.schedule( i, z_lambda( i ) ) ;
        }
      } // schedule_a_0()

      /*
      template <typename Accumulator, typename ZLambda>
      void schedule_a_1( Accumulator accumulator, ZLambda const& z_lambda ) const {
        assert( accumulator.num_columns()==num_coefs() ) ;
        for (int i=1; i<num_CORK_coefs(); ++i) {
          accumulator( glas2::all(), i ) += z_lambda( i ) ;
        }
      } // schedule_a_1()
      */
/*
      template <typename Accumulator, typename ZLambda>
      void schedule_b_0( Accumulator accumulator, ZLambda const& z_lambda ) const {
        if (num_coefs()>1)
          accumulator( glas2::all(), num_CORK_coefs() ) -= z_lambda(num_CORK_coefs()-1) ;
      } // schedule_b_0()
*/
      template <typename CoefficientMatrices, typename ZLambda>
      void schedule_b_0( CoefficientMatrices& iterator, ZLambda const& z_lambda ) const {
        if (num_coefs()>1)
          iterator.schedule( num_CORK_coefs(), -z_lambda(num_CORK_coefs()-1) ) ;
      } // schedule_b_0()

      template <typename CoefficientMatrices, typename ZLambda>
      void schedule_b_1( CoefficientMatrices& iterator, ZLambda const& z_lambda ) const {
        if (num_CORK_coefs()>1) {
          iterator.schedule( num_CORK_coefs(), -z_lambda(num_CORK_coefs()-1) ) ;
        }
      } // schedule_b_1()

/*      template <typename Accumulator, typename ZLambda>
      void schedule_b_1( Accumulator accumulator, ZLambda const& z_lambda ) const {
        if (num_CORK_coefs()>1)
          accumulator( glas2::all(), num_CORK_coefs() ) -= z_lambda(num_CORK_coefs()-1) ;
      } // schedule_b_1()
*/
    public:
      basis::monomial<I> basis_ ;
  } ; // class matrix_iterator<basis::monomial<I> >


} } // namespace CORK::matrix_iterator

#endif
