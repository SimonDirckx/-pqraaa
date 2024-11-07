//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_backend_default_backend_tensor_assign_hpp
#define glas2_backend_default_backend_tensor_assign_hpp

#include <glas2/backend/default_backend/tensor/iterator.hpp>
#include <glas2/backend/default_backend/default_backend.hpp>
#include <glas2/vector/container/static_vector.hpp>
#include <glas2/vector/container/shared_vector.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/tensor/concept/dense_tensor.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/vector/algorithm/is_equal.hpp>
#include <glas2/vector/algorithm/fill.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename To, typename From>
  typename std::enable_if< is<DenseTensor,To>::value && is<DenseTensor,From>::value, To& >::type assign( default_backend, To& to, From const& from ) {
    assert( glas2::is_equal( from.shape(), to.shape() ) );

    typedef typename To::size_type size_type ;
    size_type order = to.order() ;

    glas2::shared_vector< typename To::size_type > index( order+1 ) ;
    glas2::fill( index, 0 ) ;

    if (to.shape()(0)==0) return to ;

    for ( tensor_detail::iterator< typename To::shape_type > it(to.shape()); !it.is_end(); ++it ) {
      to( *it ) = from( *it ) ;
    }

    return to ;
  }

  template <typename To, typename From>
  typename std::enable_if< is<DenseVector,To>::value && is<DenseTensor,From>::value, To& >::type assign( default_backend, To& to, From const& from ) {
    assert( from.order() == 1 );
    assert( from.shape()(0) == to.size() );
    typedef static_vector<int,1> shape_type ;

    for (int i=0; i<to.size(); ++i) {
      to(i) = from(shape_type({i})) ;
    }

    return to ;
  }

  template <typename To, typename From>
  typename std::enable_if< is<DenseMatrix,To>::value && is<DenseTensor,From>::value, To& >::type assign( default_backend, To& to, From const& from ) {
    assert( from.order() == 2 );
    assert( from.shape()(0) == to.num_rows() );
    assert( from.shape()(1) == to.num_columns() );
    typedef static_vector<int,2> shape_type ;

    for (int i=0; i<to.num_columns(); ++i) {
      for (int j=0; j<to.num_rows(); ++j) {
        to(j,i) = from(shape_type({j,i})) ;
      }
    }

    return to ;
  }

} // namespace glas2

#endif
