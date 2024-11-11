//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_sv_aaa_my_2_hpp
#define cork_approximation_sv_aaa_my_2_hpp

#include <glas2/vector.hpp>
#include <tuple>
#include <cmath>

namespace CORK { namespace approximation {

  struct sv_aaa_2 {
    template <typename FactorMatrix>
    sv_aaa_2 apply_scaling( FactorMatrix const& ) const {
      return sv_aaa_2() ;
    }

    template <typename Combined>
    auto transform( Combined const& C ) const {
      return sv_aaa_2() ;
    }

    template <typename TestSet, typename ErrMat>
    auto operator()( TestSet const&, ErrMat const& err ) const {
       double errr=0.;
       int ind_Q=0;
       for(int i = 0;i<err.num_rows();++i){
        double nrm=glas2::norm_2(err(i,glas2::all()));
        if(errr<nrm){errr=nrm;ind_Q=i;}
       }
       return std::make_tuple( ind_Q, errr ) ;
    }
  } ; // sv_aaa_max()

} } // namespace CORK::approximation

#endif
