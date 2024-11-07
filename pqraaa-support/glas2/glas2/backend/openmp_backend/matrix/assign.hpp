//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_backend_openmp_backend_matrix_assign_hpp
#define glas2_backend_openmp_backend_matrix_assign_hpp

#include <glas2/backend/openmp_backend/openmp.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace glas2 {

  template <typename To, typename From>
  typename std::enable_if< is<DenseMatrix,To>::value && is<DenseMatrix,From>::value, To& >::type assign( openmp_backend, To& to, From const& from ) {
    assert( from.num_rows()==to.num_rows() );
    assert( from.num_columns()==to.num_columns() );

#pragma omp parallel for collapse(2)
    for (typename To::size_type i=0; i<to.num_rows(); ++i) {
      for (typename To::size_type j=0; j<to.num_columns(); ++j) {
        to(i,j) = from(i,j) ;
      }
    }
    return to ;
  }

} // namespace glas2

#endif
