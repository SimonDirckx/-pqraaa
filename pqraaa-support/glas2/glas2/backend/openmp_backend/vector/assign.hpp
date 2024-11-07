//  (C) Copyright Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_openmp_backend_vector_assign_hpp
#define glas2_openmp_backend_vector_assign_hpp

#include <glas2/backend/openmp_backend/openmp.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cassert>

namespace glas2 {

  template <typename To, typename From>
  typename std::enable_if< is<DenseVector,To>::value && is<DenseVector,From>::value, To& >::type assign( openmp, To& to, From const& from ) {
    assert( to.size() == from.size() ) ;
#pragma omp parallel for
//#pragma omp parallel for schedule(dynamic)
    for (typename To::size_type i=0; i<to.size(); ++i) { to(i) = from(i) ; }
    return to ;
  }

} // namespace glas2

#endif
