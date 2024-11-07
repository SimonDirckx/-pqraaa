//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_reval_hpp
#define cork_approximation_reval_hpp

#include <cork/approximation/reval_options.hpp>
#include <cork/approximation/reval_approximation.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <glas2/bindings.hpp>
#include <boost/numeric/bindings/lapack/driver/gesvd.hpp>
#include <boost/numeric/bindings/lapack/computational/potrf.hpp>
#include <cassert>

namespace CORK { namespace approximation {

  template <typename T, typename R, typename Z, typename F, typename W>
  T reval( T const& z, Z const& zj, F const& fj, W const& wj ) {
    assert (zj.size() == fj.size() ) ;
    assert (wj.size() == wj.size() ) ;
    T y = 0.0 ;
    T sum = 0.0 ;

    for (typename Z::size_type i=0; i<zj.size(); ++i) {
      T temp = wj(i) / (z - zj(i)) ;
      y += fj(i)*temp ;
      sum += temp ;
    }
    return y / sum ;
  } // reval()

} } // namespace CORK::approximation

#endif
