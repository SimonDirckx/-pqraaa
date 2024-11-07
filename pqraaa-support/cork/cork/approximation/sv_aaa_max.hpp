//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_sv_aaa_max_hpp
#define cork_approximation_sv_aaa_max_hpp

#include <glas2/vector.hpp>
#include <tuple>
#include <cmath>

namespace CORK { namespace approximation {

  struct sv_aaa_max {
    template <typename FactorMatrix>
    sv_aaa_max apply_scaling( FactorMatrix const& ) const {
      return sv_aaa_max() ;
    }

    template <typename Combined>
    auto transform( Combined const& C ) const {
      return sv_aaa_max() ;
    }

    template <typename TestSet, typename ErrMat>
    auto operator()( TestSet const&, ErrMat const& err ) const {
       int ind_Q = glas2::max_ind( abs( vec(err) ) ) ;
       return std::make_tuple( ind_Q % err.num_rows(), std::abs(vec(err)(ind_Q)) ) ;
    }
  } ; // sv_aaa_max()

} } // namespace CORK::approximation

#endif
