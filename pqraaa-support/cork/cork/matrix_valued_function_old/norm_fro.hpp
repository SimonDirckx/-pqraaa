//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_norm_fro_hpp
#define cork_matrix_valued_function_norm_fro_hpp

#include <cork/coefficient_matrices/norm_est.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace matrix_valued_function {


  template <typename MatrixValuedFunction>
  class norm_fro {
    public:
      typedef MatrixValuedFunction                                            matrix_valued_function_type ;
      typedef typename matrix_valued_function_type::coefficient_matrices_type coefficient_matrices_type ;
      typedef typename coefficient_matrices_type::value_type                  matrix_value_type ;
      typedef decltype( std::abs(matrix_value_type()) )                       real_type ;

    public:
      norm_fro( matrix_valued_function_type const& nep )
      : nep_( nep )
      , norm_fro_( 0 )
      {}

    public:
      void activate() {
        norm_fro_.resize( nep_.coefficient_matrices().num_matrices() ) ;

        for (typename coefficient_matrices_type::grade_type i=0; i<nep_.coefficient_matrices().num_matrices(); ++i) {
          norm_fro_( i ) = glas2::norm_fro( nep_.coefficient_matrices().sequence()[i] ) ;
        }
      }

      auto size() const {
        return nep_.basis().num_terms() ;
      }

      auto norm_coefficient( int i ) const {
        assert( norm_fro_.size()!=0 ) ;
        return norm_fro_(i) ;
      } // norm_coefficient() ;

      template <typename T>
      auto norm_for_shift( T const& shift ) const {
        assert( norm_fro_.size()!=0 ) ;
        glas2::vector< typename matrix_valued_function_type::basis_type::template value_type_for<T> > coefs( nep_.basis().num_terms() ) ;
        nep_.basis().evaluate( shift, coefs ) ;

        return inner_prod( glas2::abs(coefs), norm_fro_ ) ;
      }

      template <typename V>
      auto norm( V const& coefs ) const {
        assert( norm_fro_.size()!=0 ) ;
        return inner_prod( glas2::abs(coefs), norm_fro_ ) ;
      }

      template <typename T>
      auto upper_bound_for_shift( T const& shift ) const {
        assert( norm_fro_.size()!=0 ) ;
        return norm_for_shift( shift ) ;
      }

      template <typename V>
      auto upper_bound( V const& coefs ) const {
        assert( norm_fro_.size()!=0 ) ;
        return norm( coefs ) ;
      }

    private:
      matrix_valued_function_type const& nep_ ;
      glas2::vector< real_type >         norm_fro_ ;
  } ;

} } // namespace CORK::matrix_valued_function

#endif
