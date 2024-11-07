//  (C) Copyright Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_backend_blas_backend_blas3_assign_hpp
#define glas2_backend_blas_backend_blas3_assign_hpp

#include <glas2/backend/blas_backend/blas_backend.hpp>
#include <glas2/backend/blas_backend/blas2/is_matrix_container.hpp>
#include <glas2/matrix/expression/matrix_multiply.hpp>
#include <glas2/bindings/vector/adaptor.hpp>
#include <glas2/bindings/matrix/adaptor.hpp>
#include <boost/numeric/bindings/blas/level3/gemm.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename W, typename M, typename V>
  typename std::enable_if< is_matrix_container<W>::value && is_matrix_container<V>::value && is_matrix_container<M>::value, W >::type assign( blas_backend, W& w, matrix_multiply<M,V> const& e ) {
    assert( w.num_rows() == e.num_rows() ) ;
    assert( w.num_columns() == e.num_columns() ) ;
    boost::numeric::bindings::blas::gemm( typename W::value_type(1.0), e.matrix1(), e.matrix2(), typename W::value_type(0.0), w ) ;
    return w ;
  }

} // namespace glas2

#endif
