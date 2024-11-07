//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_sv_aaa_err_est_hpp
#define cork_matrix_valued_function_sv_aaa_err_est_hpp

#include <cork/options/aaa_stop_criterion.hpp>
#include <tuple>
#include <cmath>

namespace CORK { namespace matrix_valued_function {

  template <typename T, typename NormEst>
  class sv_aaa_err_est {
    private:
      typedef T                         value_type ;
      typedef decltype(std::abs(T()))   real_type ;

    public:
      typedef NormEst norm_est_type ;

      enum stop_criterion { MAX, WEIGHTED } ;

    public:
      sv_aaa_err_est( norm_est_type const& norm_est, glas2::range const& r, int n_fun, stop_criterion stop )
      : stop_( stop )
      , norm_est_( norm_est )
      , r_( r )
      , coefs_( norm_est_.size() )
      , factor_( n_fun )
      , norm_A_( 1.0 )
      {
   /*     switch(stop) {
          case(options::aaa_stop_criterion<real_type>::MAX):
            loew_ = NORMALIZED ;
            cvg_ = MAX ;
            break ;
          case(options::aaa_stop_criterion<real_type>::WEIGHTED):
            loew_ = WEIGHTED ;
            cvg_ = NORM ;
            break ;
        } ;
        //norm_est_.activate() ;
        */
        fill( coefs_, 0.0 ) ;
      }

    private:
      stop_criterion stop_criterion_option_2_stop( typename options::aaa_stop_criterion<real_type>::enum_type stop ) {
        stop_criterion s ;
        switch (stop) {
        case(options::aaa_stop_criterion<real_type>::MAX):
          s = MAX ;
          break;
        case(options::aaa_stop_criterion<real_type>::WEIGHTED):
          s = WEIGHTED ;
          break ;
        default:
          assert( false ) ;
        }
        return s ;
      } // stop_criterion_option_2_stop()

    public:
      sv_aaa_err_est( norm_est_type const& norm_est, glas2::range const& r, int n_fun, typename options::aaa_stop_criterion<real_type>::enum_type stop )
      : sv_aaa_err_est( norm_est, r, n_fun, stop_criterion_option_2_stop(stop) )
      {}

   
    public:
      // Scaling determines the Loewner matrix and therefore the weights.
      // MAX: divide all functions by infinity norm --> Relative error
      // WEIGHTED: scale f_i by ||A_i||
      template <typename TestSet, typename TestValues>
      void scale_functions( TestSet const& test_set, TestValues test_values ) {
        glas2::vector<real_type> f_norm( r_.size() ) ;
        glas2::vector<real_type> a_norm( r_.size() ) ;

        assert( r_.size() == test_values.num_columns() ) ;

        for (int i_fun=0; i_fun<test_values.num_columns(); ++i_fun) {
          f_norm( i_fun ) = norm_inf( test_values( glas2::all(), i_fun ) ) ;
        }

        real_type max_a_f ;
        if (stop_==WEIGHTED) {
          for (int i_mat=0; i_mat<r_.size(); ++i_mat) {
            a_norm( i_mat ) = norm_est_.norm_coefficient( r_(i_mat) ) ;
          }
          max_a_f = norm_inf( f_norm * a_norm ) ;
        std::cout << "f_norm = " << f_norm << std::endl ;
        std::cout << "a_norm = " << a_norm << std::endl ;
        }

        switch (stop_) {
          case WEIGHTED:
             norm_A_ = norm_est_.norm_for_shift( test_set(0) ) ;
             for (int i=1; i<test_set.size(); ++i) {
               norm_A_ = std::max( norm_A_, norm_est_.norm_for_shift( test_set(i) ) ) ;
             }
             break ;
          default:
            norm_A_ = 1.0 ;
        }

        real_type scal ;
        for (int i_fun=0; i_fun<test_values.num_columns(); ++i_fun) {
          switch (stop_) {
            case MAX:
              scal = norm_inf( test_values( glas2::all(), i_fun ) ) ;
              assert( scal!=0.0 ) ;
              test_values( glas2::all(), i_fun ) /= scal ;
              break;
            case WEIGHTED:
              std::cout << a_norm(glas2::slice(i_fun,a_norm.size(),test_values.num_columns())) << std::endl ;
              assert( sum(a_norm(glas2::slice(i_fun,a_norm.size(),test_values.num_columns())))!=0.0 ) ;
              //scal = sum(a_norm(glas2::slice(i_fun,a_norm.size(),test_values.num_columns()))) * max_a_f ;
              scal = a_norm(i_fun) / max_a_f ;
              test_values( glas2::all(), i_fun ) *= scal ;
              break ;
          }
          factor_(i_fun) = scal ;
        }
      } // scale_functions()

      template <typename TestValues>
      void unscale_functions( TestValues test_values ) {
        assert( test_values.num_columns()==factor_.size() ) ;

        for (int i_fun=0; i_fun<test_values.num_columns(); ++i_fun) {
          test_values( glas2::all(), i_fun ) *= factor_(i_fun) ;
        }
      } // unscale_functions()

      real_type absolute_tolerance( real_type rel_tol ) const {
        assert( stop_==MAX ) ;
        return rel_tol ;
      }

    public:
      template <typename Matrix>
      auto transform( Matrix const& map ) const {
        assert( stop_==MAX ) ;
        assert( r_.size()==map.num_rows() ) ;
        return sv_aaa_err_est( norm_est_, r_, factor_.size(), stop_ ) ;
      } // transform()

    public:
      // Return error to compare with tolerance
      // May have to be scaled.
      template <typename TestSet, typename ErrMat>
      auto operator()( TestSet const& test_set, ErrMat const& err ) const {
        assert( test_set.size() == err.num_rows() ) ;

        assert( r_.size() == err.num_columns() ) ;

        std::tuple<int,real_type> err_i ;
        fill(coefs_, 0.);

        if (stop_==MAX) {
          int ind_Q = glas2::max_ind( abs( vec(err) ) ) ;
          std::get<0>(err_i) = ind_Q % err.num_rows() ;
          std::get<1>(err_i) = std::abs(vec(err)(ind_Q)) ;
        } else {
            int i_max = 0 ;
            real_type n_max = 0.0 ;
            for (int i=0; i<err.num_rows(); ++i) {
              coefs_( r_(glas2::range(0,err.num_columns())) ) = factor_ * err( i, glas2::all() ) ;
              real_type n= norm_est_.norm( coefs_ ) ; /// norm_est_.norm_for_shift( test_set(i) ) ;
              if (n>n_max) { i_max = i ; n_max = n ; }
            }
            std::get<0>(err_i) = i_max ;
            std::get<1>(err_i) = n_max / norm_A_ ;
        }
        return err_i ;
      }

    private:
      enum stop_criterion           stop_ ;
      norm_est_type const&          norm_est_ ;
      glas2::range                  r_ ;
      glas2::vector< value_type >   coefs_ ;
      glas2::vector< real_type >    factor_ ;
      real_type                     norm_A_ ;
  } ; // sv_aaa_max()

} } // namespace CORK::matrix_valued_function

#endif
