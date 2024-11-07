//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_matrix_iterator_rational_newton_hpp
#define cork_matrix_iterator_rational_newton_hpp

#include <cork/matrix_iterator/matrix_iterator.hpp>
#include <cork/basis/rational_newton.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <cmath>
#include <limits>
#include <type_traits>

namespace CORK { namespace matrix_iterator {

  // MatrixHandler contains the matrix coefficients:
  //
  // A_i i=0,...,grade-1 contains the A-matrices
  // A_d contains B_{d-1}
  template <typename Points>
  class matrix_iterator< basis::rational_newton< Points > >
  {
    public:
      typedef basis::rational_newton< Points >  basis_type ;
      typedef typename basis_type::size_type    grade_type ;

    private:
      typedef typename basis_type::points_type  points_type ;
      typedef typename points_type::value_type  value_type ;
      typedef decltype (std::abs(value_type())) real_type ;

    public:
      matrix_iterator( basis::rational_newton< Points > const& basis )
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
        value_type factor = 1.0 ;
        if (basis_.poles()(num_CORK_coefs()-1)!=std::numeric_limits<real_type>::infinity()) factor = basis_.poles()(num_CORK_coefs()-1) ;
        for (int i=0; i<num_CORK_coefs(); ++i) {
          accumulator.schedule( i, factor * z_lambda(i) ) ;
        }
        accumulator.schedule( num_CORK_coefs(), - basis_.nodes()(num_CORK_coefs()-1) * z_lambda(num_CORK_coefs()-1) / basis_.scaling()(num_CORK_coefs()) ) ;
      }

      template <typename CoefficientMatrices, typename ZLambda>
      void schedule_a_1( CoefficientMatrices& accumulator, ZLambda const& z_lambda ) const {
        value_type factor = 1.0 ;
        if (basis_.poles()(num_CORK_coefs()-1)!=std::numeric_limits<real_type>::infinity()) factor = basis_.poles()(num_CORK_coefs()-1) ;
        for (int i=1; i<num_CORK_coefs(); ++i) {
          accumulator.schedule( i, factor * z_lambda(i) ) ;
        }
        accumulator.schedule( num_CORK_coefs(), - basis_.nodes()(num_CORK_coefs()-1) * z_lambda(num_CORK_coefs()-1) / basis_.scaling()(num_CORK_coefs()) ) ;
      }

      template <typename CoefficientMatrices, typename ZLambda>
      void schedule_b_0( CoefficientMatrices& accumulator, ZLambda const& z_lambda ) const {
        if (basis_.poles()(num_CORK_coefs()-1)!=std::numeric_limits<real_type>::infinity()) {
          for (int i=0; i<num_CORK_coefs(); ++i) {
            accumulator.schedule( i, z_lambda(i) ) ;
          }
        }
        accumulator.schedule( num_CORK_coefs(), - z_lambda(num_CORK_coefs()-1) / basis_.scaling()(num_CORK_coefs()) ) ;
      }

      template <typename CoefficientMatrices, typename ZLambda>
      void schedule_b_1( CoefficientMatrices& accumulator, ZLambda const& z_lambda ) const {
        if (basis_.poles()(num_CORK_coefs()-1)!=std::numeric_limits<real_type>::infinity()) {
          for (int i=1; i<num_CORK_coefs(); ++i) {
            accumulator.schedule( i, z_lambda(i) ) ;
          }
        }
        accumulator.schedule( num_CORK_coefs(), - z_lambda(num_CORK_coefs()-1) / basis_.scaling()(num_CORK_coefs()) ) ;
      }

    public:
      basis_type basis_ ;
  } ; // class rational_newton


} } // namespace CORK::matrix_iterator

#endif
