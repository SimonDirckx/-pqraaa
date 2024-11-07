//  (C) Copyright Karl Meerbergen 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_matrix_iterator_monomial_plus_barycentric_rational_hpp
#define cork_matrix_iterator_monomial_plus_barycentric_rational_hpp

#include <cork/utility/is_infinite.hpp>
#include <cork/matrix_iterator/matrix_iterator.hpp>
#include <cork/basis4cork/monomial_plus_barycentric_rational.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace matrix_iterator {

  // MatrixHandler contains the matrix coefficients:
  //
  // A_i i=0,...,grade-1 contains the A-matrices
  // B_i only nonzero for infinite node
  template <typename BarycentricRational>
  class matrix_iterator< basis::monomial_plus_barycentric_rational< BarycentricRational > >
  {
    public:
      typedef basis::monomial_plus_barycentric_rational< BarycentricRational >      basis_type ;
      typedef typename basis_type::size_type  grade_type ;

    public:
      matrix_iterator( basis_type const& basis )
      : basis_( basis )
      , contribute_( basis.barycentric().num_terms()-1 ) // ONLY FOR REAL CASE
      {
        assert( basis_.monomial().grade()>1 ) ;
        auto const& supports = basis.barycentric().nodes() ;
        for (int i=0; i<contribute_.size(); ++i) {
          if (glas2::imag(supports(i))==0.) contribute_[i] = true ;
          else {
            assert( i<contribute_.size()-1 ) ;
            contribute_[i] = true ;
            ++i ; contribute_[i] = false ;
          }
        }
        std::cout << "DEV: why hold a copy of basis_" << std::endl ;
      }

    public:
      grade_type num_CORK_coefs() const {
        return basis_.num_terms()-2 ;
      }

      grade_type num_coefs() const {
        return basis_.num_terms()-1 ;
      }

    private:
      grade_type num_CORK_coefs_monomial() const {
        return std::max(1,basis_.num_terms()-1) ;
      }

    public:
      template <typename Accumulator, typename ZLambda>
      void schedule_a_0( Accumulator& accumulator, ZLambda const& z_lambda ) const {
        // Schedule coefficients for degrees 1 up to d-1
        if (basis_.monomial().grade()>1) {
          // Polynomial part.
          for (int i=1; i<basis_.monomial().grade(); ++i) {
            accumulator.schedule( i, z_lambda( i-1 ) ) ;
          }

          // Constant term
          for (int i=basis_.monomial().grade()-1; i<num_CORK_coefs(); ++i) {
            if (contribute_[i-basis_.monomial().grade()+1])
              accumulator.schedule( 0, z_lambda( i ) ) ;
          }
          
          // Barycentric part
          for (int i=basis_.monomial().grade()-1; i<num_CORK_coefs(); ++i) {
            accumulator.schedule( i+2, z_lambda( i ) ) ;
          }
        } else {
          assert(false) ;
        }
      }

      template <typename Accumulator, typename ZLambda>
      void schedule_a_1( Accumulator& accumulator, ZLambda const& z_lambda ) const {
        if (basis_.monomial().grade()>1) {
          // Polynomial part.
          for (int i=2; i<basis_.monomial().grade(); ++i) {
            accumulator.schedule( i, z_lambda( i-1 ) ) ;
          }

          // Constant term
          for (int i=basis_.monomial().grade()-1; i<num_CORK_coefs(); ++i) {
            if (contribute_[i-basis_.monomial().grade()+1])
              accumulator.schedule( 0, z_lambda( i ) ) ;
          }
          
          // Barycentric part
          for (int i=basis_.monomial().grade()-1; i<num_CORK_coefs(); ++i) {
            accumulator.schedule( i+2, z_lambda( i ) ) ;
          }
        } else {
          assert(false) ;
        }
      }

      template <typename Accumulator, typename ZLambda>
      void schedule_b_0( Accumulator& accumulator, ZLambda const& z_lambda ) const {
         accumulator.schedule( basis_.monomial().grade(), -z_lambda(basis_.monomial().grade()-2) ) ;
      }

      template <typename Accumulator, typename ZLambda>
      void schedule_b_1( Accumulator& accumulator, ZLambda const& z_lambda ) const {
         if (basis_.monomial().grade()>2) accumulator.schedule( basis_.monomial().grade(), -z_lambda(basis_.monomial().grade()-2) ) ;
      }

    public:
      basis_type /*const&*/ basis_ ;
      std::vector<bool> contribute_ ;
  } ; // class monomial_plus_barycentric_rational_real


} } // namespace CORK::matrix_iterator

#endif
