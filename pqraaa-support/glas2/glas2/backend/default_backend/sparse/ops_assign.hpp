#ifndef glas2_backend_default_backend_sparse_ops_assign_hpp
#define glas2_backend_default_backend_sparse_ops_assign_hpp

#include <glas2/backend/default_backend/default_backend.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/sparse/concept/forward_coordinate_sparse_matrix.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename M, typename E>
  typename std::enable_if< is<CoordinateSparseMatrix,M>::value && is<ForwardCoordinateSparseMatrix,E>::value, M& >::type plus_assign( default_backend, M& m, E const& e ) {
    assert( m.num_rows() == e.num_rows() ) ;
    assert( m.num_columns() == e.num_columns() ) ;

    auto e_end = e.end() ;
    for ( auto it = e.begin(); it!=e_end; ++it) {
      m.push_back( it.row(), it.column(), it.value() ) ;
    }
    return m ;
  }

  template <typename To, typename From>
  typename std::enable_if< is<DenseMatrix,To>::value && is<CoordinateSparseMatrix,From>::value
                         , To&
                         >::type plus_assign( default_backend, To& to, From const& from ) {
    assert( from.num_rows()==to.num_rows() );
    assert( from.num_columns()==to.num_columns() );

    for (typename From::index_type i=0; i<from.num_nz(); ++i) {
      to( from.row(i), from.column(i) ) += from.data()(i) ;
    }
    return to ;
  }

  template <typename To, typename From>
  typename std::enable_if< is<DenseMatrix,To>::value && is<CompressedSparseMatrix,From>::value && std::is_same< typename From::orientation, row_major>::value
                         , To&
                         >::type plus_assign( default_backend, To& to, From const& from ) {
    assert( from.num_rows()==to.num_rows() );
    assert( from.num_columns()==to.num_columns() );

    for (typename From::size_type i=0; i<from.num_rows(); ++i) {
      for (typename From::index_type j=from.compressed_index()(i); j<from.compressed_index()(i+1); ++j) {
        to( i, from.direct_index()(j) ) += from.data()(j) ;
      }
    }
    return to ;
  }

} // namespace glas2

#endif
