//  (C) Copyright Karl Meerbergen 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_matrix_iterator_barycentric_rational_real_hpp
#define cork_matrix_iterator_barycentric_rational_real_hpp

#include <cork/utility/is_infinite.hpp>
#include <cork/matrix_iterator/matrix_iterator.hpp>
#include <cork/basis/barycentric_rational_real.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace matrix_iterator {

  // MatrixHandler contains the matrix coefficients:
  //
  // A_i i=0,...,grade-1 contains the A-matrices
  // B_i only nonzero for infinite node
  template <typename Weights, typename Nodes, typename I>
  class matrix_iterator< basis::barycentric_rational_real< Weights, Nodes, I > >
  {
    public:
      typedef basis::barycentric_rational_real< Weights, Nodes, I >      basis_type ;
      typedef typename basis_type::size_type  grade_type ;

    private:
      typedef typename basis_type::weights_type weights_type ;

    public:
      matrix_iterator( basis::barycentric_rational_real< Weights, Nodes, I > const& basis )
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
      template <typename Accumulator, typename ZLambda>
      void schedule_a_0( Accumulator accumulator, ZLambda const& z_lambda ) const {
        accumulator.schedule( 0, z_lambda(0) ) ;
        for (int i=1; i<num_CORK_coefs(); ++i) {
          accumulator.schedule( i, z_lambda(i) ) ;
        }
      }

      template <typename Accumulator, typename ZLambda>
      void schedule_a_1( Accumulator accumulator, ZLambda const& z_lambda ) const {
        for (int i=1; i<num_CORK_coefs(); ++i) {
          accumulator.schedule( i, z_lambda(i) ) ; //accumulator( glas2::all(), i ) += z_lambda(i) ;
        }
      }

      template <typename Accumulator, typename ZLambda>
      void schedule_b_0( Accumulator accumulator, ZLambda const& z_lambda ) const {
      }

      template <typename Accumulator, typename ZLambda>
      void schedule_b_1( Accumulator accumulator, ZLambda const& z_lambda ) const {
      }

    public:
      basis_type        basis_ ;
  } ; // class barycentric_rational_real


} } // namespace CORK::matrix_iterator

#endif
