//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_sv_aaa_norm_hpp
#define cork_matrix_valued_function_sv_aaa_norm_hpp

#include <cork/options/aaa_stop_criterion.hpp>
#include <tuple>
#include <cmath>

namespace CORK { namespace matrix_valued_function {

  template <typename T, typename FactorMatrix, typename NormEst>
  class sv_aaa_norm {
    private:
      typedef T                         value_type ;
      typedef decltype(std::abs(T()))   real_type ;

    public:
      typedef NormEst norm_est_type ;

    public:
      sv_aaa_norm( norm_est_type const& norm_est, glas2::range const& r, typename options::aaa_stop_criterion<real_type>::enum_type stop, FactorMatrix const& factor_matrix )
      : stop_( stop )
      , norm_est_( norm_est )
      , coefs_( norm_est_.size() )
      , r_( r )
      , factor_matrix_( factor_matrix )
      {
        //norm_est_.activate() ;
        fill( coefs_, 0.0 ) ;
      }

      template <typename FactorMatrixNew>
      auto set_factor_matrix( FactorMatrixNew const& factor_matrix ) const {
        return sv_aaa_norm<T, FactorMatrixNew, NormEst>( norm_est_, r_, stop_, factor_matrix ) ;
      }


    public:
      template <typename ErrMat>
      auto operator()( ErrMat const& err ) const {
         int i_max = 0 ;
         real_type n_max = 0.0 ;
         for (int i=0; i<err.num_rows(); ++i) {
           coefs_( r_ ) = multiply( factor_matrix_, err( i, glas2::all() ) ) ;
           real_type n ;
           switch (stop_) {
             case (options::aaa_stop_criterion<real_type>::WEIGHTED):
               n= norm_est_.upper_bound( coefs_ ) ;
               break;
             case (options::aaa_stop_criterion<real_type>::NORM):
               n= norm_est_.norm( coefs_ ) ;
               break;
             default:
               n = 0. ;
               assert( false ) ;
           }
           if (n>n_max) { i_max = i ; n_max = n ; }
         }

         return std::tuple( i_max, n_max ) ;
      }

    private:
      typename options::aaa_stop_criterion<real_type>::enum_type stop_ ;
      norm_est_type const&                                       norm_est_ ;
      glas2::vector< value_type >                                coefs_ ;
      glas2::range                                               r_ ;
      FactorMatrix                                               factor_matrix_ ;
  } ; // sv_aaa_max()

} } // namespace CORK::matrix_valued_function

#endif
