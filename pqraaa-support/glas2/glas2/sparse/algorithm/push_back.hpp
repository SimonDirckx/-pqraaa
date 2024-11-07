#ifndef glas2_sparse_algorithm_push_back_hpp
#define glas2_sparse_algorithm_push_back_hpp

#include <glas2/concept/is.hpp>
#include <glas2/sparse/concept/coordinate_sparse_matrix.hpp>
#include <glas2/sparse/concept/compressed_sparse_matrix.hpp>
#include <glas2/sparse/view/sparse_upper_view.hpp>
#include <glas2/sparse/concept/sparse_matrix.hpp>
#include <type_traits>

namespace glas2 {

  namespace detail {

    struct no_filter {
      template <typename COO, typename I, typename C, typename T>
      void operator() ( COO& coo, I row, C col, T const& v ) const {
        coo.push_back( row, col, v ) ;
      }
    } ;

    struct upper_filter {
      template <typename COO, typename I, typename C, typename T>
      void operator() ( COO& coo, I row, C col, T const& v ) const {
        if (row<=col) coo.push_back( row, col, v ) ;
      }
    } ;

  } // namespace detail

  template <typename C, typename M, typename R1, typename R2, typename Filter=detail::no_filter>
  C& push_back( C& coo, M const& input, R1 const& rows, R2 const& columns, Filter const& filter=detail::no_filter() ) {
    for (typename M::size_type i=0; i<input.num_nz(); ++i) {
      filter( coo,  rows(input.row(i)), columns(input.column(i)), input.data()(i) ) ;
      //coo.push_back( rows(input.row(i)), columns(input.column(i)), input.data()(i) ) ;
    }
    return coo ;
  }

  // Coordinate
  template <typename C, typename M, typename Filter=detail::no_filter>
  typename std::enable_if< is<CoordinateSparseMatrix, M>::value, C&>::type push_back( C& coo, M const& input, Filter const& filter=detail::no_filter() ) {
    for (typename M::size_type i=0; i<input.num_nz(); ++i) {
      filter( coo, input.row(i), input.column(i), input.data()(i) ) ;
      //coo.push_back( input.row(i), input.column(i), input.data()(i) ) ;
    }
    return coo ;
  }

  template <typename C, typename M>
  C& push_back( C& coo, sparse_upper_view<M> const& input ) {
    return push_back( coo, input.matrix(), detail::upper_filter() ) ;
  }

  // Compressed
  namespace detail {
    template <typename C, typename M, typename Filter>
    void push_back_crs( column_major , C& coo, M const& input, Filter const& filter ) {
      for (typename M::size_type i=0; i<input.num_columns(); ++i) {
        for (typename M::size_type j=input.compressed_index()(i); j<input.compressed_index()(i+1); ++j) {
          filter( coo, input.direct_index()(j), i, input.data()(j) ) ;
          //coo.push_back( input.direct_index()(j), i, input.data()(j) ) ;
        }
      }
    }

    template <typename C, typename M, typename Filter>
    void push_back_crs( row_major , C& coo, M const& input, Filter const& filter ) {
      for (typename M::size_type i=0; i<input.num_rows(); ++i) {
        for (typename M::size_type j=input.compressed_index()(i); j<input.compressed_index()(i+1); ++j) {
          filter( coo, i, input.direct_index()(j), input.data()(j) ) ;
          //coo.push_back( i, input.direct_index()(j), input.data()(j) ) ;
        }
      }
    }
  } // namespace detail

  template <typename C, typename M, typename Filter=detail::no_filter>
  typename std::enable_if< is<CompressedSparseMatrix, M>::value, C&>::type push_back( C& coo, M const& input, Filter const& filter=detail::no_filter() ) {
    detail::push_back_crs( typename M::orientation(), coo, input, filter ) ;
    return coo ;
  }

}

#endif
