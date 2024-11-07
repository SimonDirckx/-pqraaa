//  (C) Copyright Karl Meerbergen 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_matrix_iterator_barycentric_rational_real_strong_hpp
#define cork_matrix_iterator_barycentric_rational_real_strong_hpp

#include <cork/utility/is_infinite.hpp>
#include <cork/matrix_iterator/matrix_iterator.hpp>
#include <cork/basis/barycentric_rational_real_strong.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace matrix_iterator {

  // MatrixHandler contains the matrix coefficients:
  //
  // Terms have the form: A_i + z * B_i
  // A_i i=0,...,grade-1 contains the A-matrices: takes numbers 1 ... num_CORK_coefs() (A_0 is not used)
  // B_i i=0,...,grade-1 contains the B-matrices: takes numbers num_CORK_coefs()+1 ... 2*num_CORK_coefs()
  // No infinite nodes allowed.
  template <typename Weights, typename Nodes, typename I>
  class matrix_iterator< basis::barycentric_rational_real_strong< Weights, Nodes, I > >
  {
    public:
      typedef basis::barycentric_rational_real_strong< Weights, Nodes, I >      basis_type ;
      typedef typename basis_type::size_type  grade_type ;

    private:
      typedef typename basis_type::weights_type weights_type ;

    public:
      matrix_iterator( basis::barycentric_rational_real_strong< Weights, Nodes, I > const& basis )
      : basis_( basis )
      , real_basis_( std::imag(basis_.nodes()(basis_.nodes().size()-1))==0.0 )
      {}

    public:
      grade_type num_CORK_coefs() const {
        return basis_.num_terms()-1 ;
      }

      grade_type num_coefs() const {
        return basis_.num_terms() ; // for A_0...A_d, and B_1...B_d
      }

    public:
      template <typename Accumulator, typename ZLambda>
      void schedule_a_0( Accumulator accumulator, ZLambda const& z_lambda ) const {
        schedule_a_1( accumulator, z_lambda ) ;
        // Last but first basis function may have a contribution to term 0.
        if (real_basis_)
          accumulator.schedule( num_CORK_coefs(), z_lambda(0) ) ;
        else
          accumulator.schedule( num_CORK_coefs()-1, z_lambda(0) ) ;
      }

      template <typename Accumulator, typename ZLambda>
      void schedule_a_1( Accumulator accumulator, ZLambda const& z_lambda ) const {
        for (int i=1; i<num_CORK_coefs()-1; ++i) {
          assert( !is_infinite(basis_.nodes()(i-1)) ) ;
          accumulator.schedule( i, z_lambda(i) ) ; //accumulator( glas2::all(), i ) += z_lambda(i) ;
        }
        // Last basis functions
        if (real_basis_) {
          accumulator.schedule( num_CORK_coefs()-1, z_lambda(num_CORK_coefs()-1) ) ;
          for (int i=1; i<num_CORK_coefs(); ++i) {
            accumulator.schedule( num_CORK_coefs(), -z_lambda(i) ) ;
          }
        } else {
          accumulator.schedule( num_CORK_coefs(), z_lambda(num_CORK_coefs()-1) ) ; //accumulator( glas2::all(), i ) += z_lambda(i) ;
          for (int i=0; i<basis_.nodes().size()-2; ++i) {
            accumulator.schedule( num_CORK_coefs()-1, -z_lambda(i+1) ) ;
            if (glas2::imag(basis_.nodes()(i))!=0.0) ++i ;
          }
        }
      }

      template <typename Accumulator, typename ZLambda>
      void schedule_b_0( Accumulator accumulator, ZLambda const& z_lambda ) const {
        schedule_b_1( accumulator, z_lambda ) ;
        // Last but first basis function may have a contribution to term 0.
 //       accumulator.schedule( 2*num_CORK_coefs()-1, -z_lambda(0) ) ;
      }

      template <typename Accumulator, typename ZLambda>
      void schedule_b_1( Accumulator accumulator, ZLambda const& z_lambda ) const {
 /*       for (int i=1; i<num_CORK_coefs()-1; ++i) {
          assert( !is_infinite(basis_.nodes()(i-1)) ) ;
          accumulator.schedule( num_CORK_coefs()+i, -z_lambda(i) ) ; //accumulator( glas2::all(), i ) += z_lambda(i) ;
        }
        accumulator.schedule( num_CORK_coefs()*2, -z_lambda(num_CORK_coefs()-1) ) ; //accumulator( glas2::all(), i ) += z_lambda(i) ;
        for (int i=0; i<basis_.nodes().size()-2; ++i) {
          accumulator.schedule( 2*num_CORK_coefs()-1, z_lambda(i+1) ) ;
          if (glas2::imag(basis_.nodes()(i))!=0.0) ++i ;
        }*/
      }

    public:
      basis_type basis_ ;
      bool       real_basis_ ;
  } ; // class barycentric_rational_real_strong


} } // namespace CORK::matrix_iterator

#endif
