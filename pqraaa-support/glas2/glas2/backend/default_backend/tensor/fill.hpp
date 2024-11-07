//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_backend_default_backend_tensor_algorithm_fill_hpp
#define glas2_backend_default_backend_tensor_algorithm_fill_hpp

#include <glas2/backend/default_backend/default_backend.hpp>
#include <glas2/tensor/concept/contiguous_dense_tensor.hpp>
#include <glas2/tensor/concept/strided_dense_tensor.hpp>
#include <glas2/vector/container/shared_vector.hpp>
#include <glas2/concept/is.hpp>
#include <algorithm>
#include <type_traits>
#include <vector>

namespace glas2 {

  template <typename X, typename T>
  typename std::enable_if< is<ContiguousDenseTensor,X>::value, X >::type fill( current_backend, X x, T const& value ) {
    typedef typename X::size_type size_type ;
    //std::vector<size_type> index( x.order()+1 ) ;
    //std::fill( index.begin(), index.end(), 0 ) ;

    size_type nnz = prod( x.shape() ) ;
    for (size_type i=0; i<nnz; ++i) {
      x[i] = value ;
      //index[0]++ ; int k=0; while ( k<x.order() && index[k]==x.shape()(k) ) { index[k] = 0; ++k; index[k]++ ; }
    }

    return x ;
  }

  template <typename X, typename T>
  typename std::enable_if< is<StridedDenseTensor,X>::value, X >::type fill( current_backend, X x, T const& value ) {
    typedef typename X::size_type size_type ;
    glas2::shared_vector<size_type> index( x.order()+1 ) ;
    auto order = x.order() ;
    auto index_r = index( range(0,order) ) ;

    fill( index, 0 ) ;
    while (index(order)==0) {
      x( index_r ) = value ;     
      ++index(0) ;
      int j=0;
      bool done = index(j)<x.shape()(j) ;
      while (!done) {
        index(j) = 0 ;
        ++ j ;
        ++ index(j) ;
        if (j<order) done = index(j)<x.shape()(j) ;
        else done = true ;
      }
    }

    return x ;
  }

} // namespace glas2

#endif
