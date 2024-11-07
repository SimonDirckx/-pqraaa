#ifndef glas2_sparse_algorithm_transpose_hpp
#define glas2_sparse_algorithm_transpose_hpp

#include <glas2/concept/is.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/sparse/type/crs_structure.hpp>
#include <glas2/sparse/type/coo_structure.hpp>
#include <type_traits>

namespace glas2 {

  template <typename M>
  typename std::enable_if< is<CompressedSparseMatrix,M>::value && std::is_same< typename M::orientation, row_major >::value
                         , crs_structure<typename M::value_type, column_major, typename M::index_type, typename M::size_type, M::index_base >
                         >::type transpose( M m ) {
    return crs_structure<typename M::value_type, column_major, typename M::index_type, typename M::size_type, M::index_base>( m.num_columns(), m.num_rows(), m.num_nz(), m.data(), m.compressed_index(), m.direct_index() ) ;
  }

  template <typename M>
  typename std::enable_if< is<CompressedSparseMatrix,M>::value && std::is_same< typename M::orientation, column_major >::value
                         , crs_structure<typename M::value_type, row_major, typename M::index_type, typename M::size_type, M::index_base >
                         >::type transpose( M m ) {
    return crs_structure<typename M::value_type, row_major, typename M::index_type, typename M::size_type, M::index_base>( m.num_columns(), m.num_rows(), m.num_nz(), m.data(), m.compressed_index(), m.direct_index() ) ;
  }

  template <typename M>
  typename std::enable_if< is<CoordinateSparseMatrix,M>::value
                         , coo_structure<typename M::value_type, typename M::index_type, typename M::size_type, M::index_base >
                         >::type transpose( M m ) {
    return coo_structure<typename M::value_type, typename M::index_type, typename M::size_type, M::index_base>( m.num_columns(), m.num_rows(), m.num_nz(), m.storage_data(), m.storage_columns(), m.storage_rows() ) ;
  }

} // namespace glas2

#endif
