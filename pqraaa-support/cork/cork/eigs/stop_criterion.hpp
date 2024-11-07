//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_eigs_stop_criterion_hpp
#define cork_eigs_stop_criterion_hpp

#include <cork/options/stop_criterion.hpp>
#include <cork/matrix_valued_function/norm_est.hpp>
#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <iosfwd>

namespace CORK { namespace eigs {

  template <typename T, typename NEP, typename StopCriterion>
  class stop_criterion {
    private:
      typedef T                                 value_type ;
      typedef typename StopCriterion::real_type real_type ;

    public:
      stop_criterion( NEP const& nep, matrix_valued_function::norm_est<NEP>& norm_est, StopCriterion const& stop )
      : nep_( nep )
      , stop_criterion_( stop )
      , resid_( 0 )
      , norm_est_( norm_est )
      {
        if (stop_criterion_.choice()==StopCriterion::NEP || stop_criterion_.choice()==StopCriterion::MIXED) {
          norm_est_.activate() ;
          resid_.resize( nep_.size() ) ;
        }
      }
      
    public:
      bool eigenvector_needed() const {
        return resid_.size()!=0 ;
      }

      template <typename Vec>
      bool test( value_type const& value, Vec const& eig_vec, real_type const& rks_resid, real_type& nep_resid, real_type const& norm_h ) const {
        real_type tolerance ;
        switch (stop_criterion_.choice()) {
          case StopCriterion::KRYLOV:
            tolerance = stop_criterion_.absolute_tolerance() + norm_h * std::max( std::numeric_limits<real_type>::epsilon(), stop_criterion_.relative_tolerance() ) ;
            return rks_resid <= tolerance ;
          case StopCriterion::MIXED:
            tolerance = stop_criterion_.absolute_tolerance() + norm_h * std::max( std::numeric_limits<real_type>::epsilon(), stop_criterion_.relative_tolerance() ) ;
            if (rks_resid > tolerance) return false ;
          case StopCriterion::NEP:
            fill( resid_, 0.0 ) ;
            nep_.multiply_add( value, eig_vec, resid_ ) ;
            nep_resid = norm_2(resid_) ;
            tolerance = stop_criterion_.absolute_tolerance() + norm_est_.norm_for_shift(value) * std::max( std::numeric_limits<real_type>::epsilon(), stop_criterion_.relative_tolerance() ) ;
            return nep_resid <= tolerance ;
          case StopCriterion::INVARIANT_PAIR:
            assert(false) ;
            return false ;
        }
      } // test()

    private:
      NEP const&                             nep_ ;
      StopCriterion                          stop_criterion_ ;
      glas2::vector< value_type >            resid_ ;
      matrix_valued_function::norm_est<NEP>& norm_est_ ;
  } ; // class stop_criterion

} } // namespace CORK::eigs


#endif
