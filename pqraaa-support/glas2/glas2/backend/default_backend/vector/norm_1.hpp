//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_backend_default_backend_vector_algorithm_norm_1_hpp
#define glas2_backend_default_backend_vector_algorithm_norm_1_hpp

#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/backend/default_backend/default_backend.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cmath>
#ifdef GLAS_OPENMP
//#include <omp.h>
#endif

namespace glas2 {

  template <typename X>
  typename std::enable_if< is<DenseVector,X>::value
                         , decltype( std::abs(typename X::value_type()) )
                         >::type norm_1( default_backend, X const& x ) {
    decltype( std::abs(typename X::value_type()) ) sum = 0 ;
#ifdef GLAS_OPENMP
#pragma omp parallel for reduction (+:sum)
#endif
    for (typename X::size_type i=0; i<x.size(); ++i) {
//      std::cout << omp_get_thread_num() << " / " << omp_get_num_threads() << std::endl ;
      sum += std::abs( x(i) ) ;
    }
    return sum ;
  }
} // namespace glas2

#endif
