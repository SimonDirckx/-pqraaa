//  (C) Copyright Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_blas_openmp_backend_vector_assign_hpp
#define glas2_blas_openmp_backend_vector_assign_hpp

#include <glas2/backend/blas_backend/blas1.hpp>
#include <glas2/backend/blas_backend/blas2.hpp>
#include <glas2/backend/blas_openmp_backend/blas_openmp.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cassert>

namespace glas2 {

  template <typename To, typename From>
  typename std::enable_if< is<DenseVector,To>::value && is<DenseVector,From>::value, To& >::type assign( blas_openmp, To& to, From const& from ) {
    assert( to.size() == from.size() ) ;
#ifdef _OPENMP
    int n_threads = omp_get_num_threads() ;
#else
    int n_threads = 1 ;
#endif
    typename To::size_type step = (to.size()+n_threads-1) / n_threads ;
#pragma omp parallel for
    for (int i=0; i<n_threads; ++i) {
      glas2::range r( step*i, std::min(to.size(),step*i+step) ) ;
      auto to_r = to( r ) ;
      assign( blas_backend(), to_r, from(r) ) ;
    }
    return to ;
  }

} // namespace glas2

#endif
