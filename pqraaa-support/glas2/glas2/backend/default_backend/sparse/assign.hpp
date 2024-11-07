//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_backend_default_backend_sparse_assign_hpp
#define glas2_backend_default_backend_sparse_assign_hpp

#include <glas2/backend/default_backend/default_backend.hpp>
#include <glas2/backend/default_backend/vector/fill.hpp>
#include <glas2/backend/default_backend/matrix/fill.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/sparse/concept/forward_coordinate_sparse_matrix.hpp>
#include <glas2/sparse/concept/compressed_sparse_matrix.hpp>
#include <glas2/sparse/concept/sparse_vector.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename To, typename From>
  typename std::enable_if< is<DenseMatrix,To>::value && is<CoordinateSparseMatrix,From>::value
                         , To&
                         >::type assign( default_backend, To& to, From const& from ) {
    assert( from.num_rows()==to.num_rows() );
    assert( from.num_columns()==to.num_columns() );

    fill( default_backend(), to, 0.0 ) ;

    for (typename From::index_type i=0; i<from.num_nz(); ++i) {
      to( from.row(i), from.column(i) ) += from.data()(i) ;
    }
    return to ;
  }

  template <typename To, typename From>
  typename std::enable_if< is<DenseVector,To>::value && is<SparseVector,From>::value
                         , To&
                         >::type assign( default_backend, To& to, From const& from ) {
    assert( from.size()==to.size() );

    fill( default_backend(), to, 0.0 ) ;

    for (typename From::size_type i=0; i<from.num_nz(); ++i) {
      to( from.index()(i) ) += from.data()(i) ;
    }
    return to ;
  }

  template <typename To, typename From>
  typename std::enable_if< is<DenseMatrix,To>::value && is<CompressedSparseMatrix,From>::value && std::is_same< typename From::orientation, row_major>::value
                         , To&
                         >::type assign( default_backend, To& to, From const& from ) {
    assert( from.num_rows()==to.num_rows() );
    assert( from.num_columns()==to.num_columns() );

    fill( default_backend(), to, 0.0 ) ;

    for (typename From::size_type i=0; i<from.num_rows(); ++i) {
      for (typename From::index_type j=from.compressed_index()(i); j<from.compressed_index()(i+1); ++j) {
        to( i, from.direct_index()(j) ) += from.data()(j) ;
      }
    }
    return to ;
  }

} // namespace glas2

#endif
