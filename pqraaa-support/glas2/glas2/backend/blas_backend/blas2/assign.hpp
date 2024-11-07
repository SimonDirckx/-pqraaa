//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_backend_blas_backend_blas2_assign_hpp
#define glas2_backend_blas_backend_blas2_assign_hpp

#include <glas2/backend/blas_backend/blas_backend.hpp>
#include <glas2/backend/blas_backend/blas1/is_vector_container.hpp>
#include <glas2/backend/blas_backend/blas2/is_matrix_container.hpp>
#include <glas2/matrix/expression/matrix_vector_multiply.hpp>
#include <glas2/bindings/vector/adaptor.hpp>
#include <glas2/bindings/matrix/adaptor.hpp>
#include <glas2/bindings/matrix/expression.hpp>
#include <boost/numeric/bindings/blas/level2/gemv.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename W, typename M, typename V>
  typename std::enable_if< is_vector_container<W>::value && is_vector_container<V>::value && is_matrix_container<M>::value, W >::type assign( blas_backend, W& w, matrix_vector_multiply<M,V> const& e ) {
    assert( w.size() == e.size() ) ;
    boost::numeric::bindings::blas::gemv( typename W::value_type(1.0), e.matrix(), e.vector(), typename W::value_type(0.0), w ) ;
    return w ;
  }

} // namespace glas2

#endif
