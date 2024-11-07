//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_norm_est_hpp
#define cork_matrix_valued_function_norm_est_hpp

#include <cork/coefficient_matrices/norm_est.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>
#include <memory>

namespace CORK { namespace matrix_valued_function {


  template <typename MatrixValuedFunction>
  class norm_est {
    public:
      typedef MatrixValuedFunction                                                        matrix_valued_function_type ;
      typedef typename matrix_valued_function_type::coefficient_matrices_type::value_type value_type ;

    protected:
      typedef coefficient_matrices::norm_est< value_type >                                norm_est_type ;

    public:
      norm_est( matrix_valued_function_type const& nep )
      : nep_( nep )
      {}

    public:
      void activate() {
#ifndef NDEBUG
        if (norm_est_.get()) std::cerr << "WARNING: CORK::matrix_valued_function::norm_est::activate(): is already activated" << std::endl ;
#endif
        if (!norm_est_.get() ) norm_est_.reset( new norm_est_type( nep_.coefficient_matrices() ) ) ;
      }

    public:
      auto const& nep() const { return nep_ ; }

      auto size() const {
        return nep_.basis().num_terms() ;
      }

      auto const& coefficients_norm_est() const { return *norm_est_ ; }

    public:
      auto norm_coefficient( int i ) const {
        assert( norm_est_.get() ) ;
        return norm_est_->norm_coefficient(i) ;
      } // norm_coefficient() ;

      template <typename T>
      auto norm_for_shift( T const& shift ) const {
        assert( norm_est_.get() ) ;
        glas2::vector< typename matrix_valued_function_type::basis_type::template value_type_for<T> > coefs( nep_.basis().num_terms() ) ;
        nep_.basis().evaluate( shift, coefs ) ;

        return norm_est_->norm( coefs ) ;
      }

      template <typename V>
      auto norm( V const& coefs ) const {
        assert( norm_est_.get() ) ;
        static_assert( glas2::is< glas2::Vector, V >::value ) ;
        return norm_est_->norm( coefs ) ;
      }

      template <typename T>
      auto upper_bound_for_shift( T const& shift ) const {
        assert( norm_est_.get() ) ;
        glas2::vector< typename matrix_valued_function_type::basis_type::template value_type_for<T> > coefs( nep_.basis().num_terms() ) ;
        nep_.basis().evaluate( shift, coefs ) ;

        return norm_est_->upper_bound( coefs ) ;
      }

      template <typename V>
      auto upper_bound( V const& coefs ) const {
        assert( norm_est_.get() ) ;
        static_assert( glas2::is< glas2::Vector, V >::value ) ;
        return norm_est_->upper_bound( coefs ) ;
      }

    private:
      matrix_valued_function_type const& nep_ ;
      std::shared_ptr<norm_est_type>     norm_est_ ;
  } ;

} } // namespace CORK::matrix_valued_function

#endif
