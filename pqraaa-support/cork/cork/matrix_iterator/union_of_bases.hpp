//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_matrix_polynomial_union_of_bases_matrix_polynomial_hpp
#define cork_matrix_polynomial_union_of_bases_matrix_polynomial_hpp

#include <cork/coefficient_matrices/selection.hpp>
#include <cork/matrix_iterator/matrix_iterator.hpp>
#include <cork/basis/union_of_bases.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace matrix_iterator {

  // MatrixHandler contains the matrix coefficients:
  //
  // A_i i=0,...,grade-1 contains the A-matrices
  // A_d contains B_{d-1}
  template <typename B1, typename B2>
  class matrix_iterator< basis::union_of_bases<B1,B2> >
  {
    public:
      typedef basis::union_of_bases< B1, B2 >         basis_type ;
      typedef typename basis_type::size_type  grade_type ;

    public:
      matrix_iterator( basis::union_of_bases< B1, B2 > const& basis )
      : basis_( basis )
      , iterator_1_( basis.basis_1() )
      , iterator_2_( basis.basis_2() )
      , range_1_( 0, iterator_1_.num_coefs() )
      , range_2_( iterator_2_.num_coefs() )
      , range_2_CORK_( iterator_2_.num_CORK_coefs() )
      {
        range_2_(0) = 0 ;
        range_2_( glas2::range_from_end(1,0) ) = glas2::range(iterator_1_.num_coefs(),num_coefs()) ;
        range_2_CORK_(0) = 0 ;
        range_2_CORK_( glas2::range_from_end(1,0) ) = glas2::range(iterator_1_.num_CORK_coefs(),num_CORK_coefs()) ;
      }

    public:
      grade_type num_CORK_coefs() const {
        return iterator_1_.num_CORK_coefs() + iterator_2_.num_CORK_coefs() - 1 ;
      }

      grade_type num_coefs() const {
        return iterator_1_.num_coefs() + iterator_2_.num_coefs() - 1 ;
      }

    public:
      template <typename CoefficientMatrices, typename ZLambda>
      void schedule_a_0( CoefficientMatrices& iterator, ZLambda const& z_lambda ) const {
        auto it1 = iterator.range( range_1_ ) ;
        iterator_1_.schedule_a_0( it1, z_lambda ) ;
        typename decltype(range_2_CORK_)::base_type range_2( range_2_CORK_ ) ;
        auto it2 = iterator.range(range_2_) ;
        iterator_2_.schedule_a_0( it2, [&range_2,&z_lambda](auto i) { return z_lambda(range_2(i)) ;} ) ;
      }

      // Different flavours of scheduling, with different ranges
      // The first block is not taken, so, no need to exclude A_0
      template <typename CoefficientMatrices, typename ZLambda>
      void schedule_a_1( CoefficientMatrices& iterator, ZLambda const& z_lambda ) const {
        auto it1 = iterator.range(range_1_) ;
        iterator_1_.schedule_a_1( it1, z_lambda ) ;
        typename decltype(range_2_CORK_)::base_type range_2( range_2_CORK_ ) ;
        auto it2 = iterator.range(range_2_) ;
        iterator_2_.schedule_a_1( it2, [&range_2,&z_lambda](auto i) { return z_lambda(range_2(i)) ;} ) ;
      }

      template <typename CoefficientMatrices, typename ZLambda>
      void schedule_b_0( CoefficientMatrices& iterator, ZLambda const& z_lambda ) const {
        auto it1 = iterator.range(range_1_) ;
        iterator_1_.schedule_b_0( it1, z_lambda ) ;
        typename decltype(range_2_CORK_)::base_type range_2( range_2_CORK_ ) ;
//        z_lambda_type_0<ZLambda> z_lambda_0( z_lambda ) ;
        // No need to take block 0 again:
        auto it2 = iterator.range(range_2_) ;
        iterator_2_.schedule_b_0( it2, [&range_2,&z_lambda](auto i) { return z_lambda(range_2(i)) ;} ) ;
      }

      template <typename CoefficientMatrices, typename ZLambda>
      void schedule_b_1( CoefficientMatrices& iterator, ZLambda const& z_lambda ) const {
        auto it1 = iterator.range(range_1_) ;
        iterator_1_.schedule_b_1( it1, z_lambda ) ;
        typename decltype(range_2_CORK_)::base_type range_2( range_2_CORK_ ) ;
        auto it2 = iterator.range(range_2_) ;
        iterator_2_.schedule_b_1( it2, [&range_2,&z_lambda](auto i) { return z_lambda(range_2(i)) ;} ) ;
      }

    public:
      basis_type                                           basis_ ;
      matrix_iterator< typename basis_type::basis_1_type > iterator_1_ ;
      matrix_iterator< typename basis_type::basis_2_type > iterator_2_ ;
      glas2::range                                         range_1_ ;
      glas2::shared_vector<int>                            range_2_ ;
      glas2::shared_vector<int>                            range_2_CORK_ ;
  } ; // class union_of_bases


} } // namespace CORK::matrix_iterator

#endif
