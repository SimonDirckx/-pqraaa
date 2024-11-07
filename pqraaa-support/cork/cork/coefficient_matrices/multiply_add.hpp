//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_coefficient_matrices_multiply_add_hpp
#define cork_coefficient_matrices_multiply_add_hpp

#include <cassert>
#include <type_traits>

namespace CORK { namespace coefficient_matrices {

  template <typename CoefficientMatrices>
  struct multiply_add {
    template <typename Coefficients, typename Range>
    static void apply( CoefficientMatrices const& coefficient_matrices, Coefficients const& coefs, Range const& range, X const& x, W w ) {
      for (typename Range::size_type i=0; i<range.size(); ++i) {
        assert( range(i)<coefficient_matrices.num_matrices() && range(i)>=0 ) ;
        coefficient_matrices.multiply_add( range(i), x, w ) ;
      }
    }
  } ; // multiply_add


} } // namespace CORK::coefficient_matrices

#endif
