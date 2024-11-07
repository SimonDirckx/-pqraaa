//  (C) Copyright Karl Meerbergen 2022.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_sv_aaa_err_est_base_hpp
#define cork_matrix_valued_function_sv_aaa_err_est_base_hpp

#include <cork/options/aaa_stop_criterion.hpp>
#include <tuple>
#include <cmath>

namespace CORK { namespace matrix_valued_function {

  template <typename T, typename NormEst, typename Factor>
  class sv_aaa_err_est_base {
    private:
      typedef T                                                          value_type ;
      typedef decltype(std::abs(T()))                                    real_type ;
      typedef typename options::aaa_stop_criterion<real_type>::enum_type stop_type ;

    public:
      typedef NormEst norm_est_type ;

      enum stop_criterion { MAX, WEIGHTED } ;

    public:
      sv_aaa_err_est_base( norm_est_type const& norm_est, glas2::range const& r, Factor const& factor, stop_type stop )
      : stop_( stop )
      , norm_est_( norm_est )
      , r_( r )
      , coefs_( r.size() )
      , factor_( factor )
      {}

    public:
      sv_aaa_err_est( norm_est_type const& norm_est, glas2::range const& r, int n_fun, stop_type stop )
      : sv_aaa_err_est( norm_est, r, n_fun, stop_criterion_option_2_stop(stop) )
      {}

    public:
      // Return error to compare with tolerance
      // May have to be scaled.
      template <typename TestSet, typename ErrMat>
      auto operator()( TestSet const& test_set, ErrMat const& err ) const {
        assert( test_set.size() == err.num_rows() ) ;

        assert( r_.size() == err.num_columns() ) ;

        std::tuple<int,real_type> err_i ;
        fill(coefs_, 0.);

        switch (stop_) {
        case (stop_type::MAX):
          int ind_Q = glas2::max_ind( abs( vec(err) ) ) ;
          std::get<0>(err_i) = ind_Q % err.num_rows() ;
          std::get<1>(err_i) = std::abs(vec(err)(ind_Q)) ;
          break ;
        case (stop_type::WEIGHTED):
          int i_max = 0 ;
          real_type n_max = 0.0 ;
          for (int i=0; i<err.num_rows(); ++i) {
            coefs_( r_(glas2::range(0,err.num_columns())) ) =  multiply(factor_, err( i, glas2::all() ) ) ;
            real_type n= norm_est_.norm( coefs_ ) ; /// norm_est_.norm_for_shift( test_set(i) ) ;
            if (n>n_max) { i_max = i ; n_max = n ; }
          }
          std::get<0>(err_i) = i_max ;
          std::get<1>(err_i) = n_max / norm_A_ ;
          break ;

        case default:
          assert(false) ;
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
  } ; // class sv_aaa_err_est_bas)

} } // namespace CORK::matrix_valued_function

#endif
