//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_norm_1_hpp
#define cork_matrix_valued_function_norm_1_hpp

#include <cork/coefficient_matrices/norm_est.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace matrix_valued_function {


  template <typename MatrixValuedFunction>
  class norm_1 {
    public:
      typedef MatrixValuedFunction matrix_valued_function_type ;

    public:
      norm_1( matrix_valued_function_type const& nep )
      : nep_( nep )
      , norm_1_( nep.coefficient_matrices().num_matrices() )
      {
      }

      ~norm_1() {
        delete norm_1_ ;
      }

    public:
      void activate() {
        if (!norm_1_) norm_1_ = new norm_est_type( nep_.coefficient_matrices() ) ;
      }

      auto size() const {
        return nep_.basis().num_terms() ;
      }

      auto norm_coefficient( int i ) const {
        assert( norm_1_ ) ;
        return norm_1_->norm_coefficient(i) ;
      } // norm_coefficient() ;

      template <typename T>
      auto norm_for_shift( T const& shift ) const {
        assert( norm_est_ ) ;
        glas2::vector< typename matrix_valued_function_type::basis_type::template value_type_for<T> > coefs( nep_.basis().num_terms() ) ;
        nep_.basis().evaluate( shift, coefs ) ;

        return norm_est_->norm( coefs ) ;
      }

      template <typename V>
      auto norm( V const& coefs ) const {
        assert( norm_est_ ) ;
        return norm_est_->norm( coefs ) ;
      }

      template <typename T>
      auto upper_bound_for_shift( T const& shift ) const {
        assert( norm_est_ ) ;
        glas2::vector< typename matrix_valued_function_type::basis_type::template value_type_for<T> > coefs( nep_.basis().num_terms() ) ;
        nep_.basis().evaluate( shift, coefs ) ;

        return norm_est_->upper_bound( coefs ) ;
      }

      template <typename V>
      auto upper_bound( V const& coefs ) const {
        assert( norm_est_ ) ;
        return norm_est_->upper_bound( coefs ) ;
      }

    private:
      matrix_valued_function_type const& nep_ ;
      norm_est_type*                     norm_est_ ;
  } ;

} } // namespace CORK::matrix_valued_function

#endif
