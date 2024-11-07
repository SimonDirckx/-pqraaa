//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_sparse_container_crs_hpp
#define glas2_sparse_container_crs_hpp

#include <glas2/vector/container/shared_vector.hpp>
#include <glas2/sparse/algorithm/assign.hpp>
#include <glas2/sparse/algorithm/is_sorted.hpp>
#include <glas2/sparse/concept/compressed_sparse_matrix.hpp>
#include <glas2/sparse/concept/coordinate_sparse_matrix.hpp>
#include <glas2/sparse/type/crs_structure.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

/*  namespace detail {

    template <typename Orientation>
    struct size_of_compressed_index
    {} ;

    template <>
    struct size_of_compressed_index< ::glas2::row_major >
    {
      template <typename S>
      static S apply( S n, S m ) { return n+1 ; }
    } ;

    template <>
    struct size_of_compressed_index< ::glas2::column_major >
    {
      template <typename S>
      static S apply( S n, S m ) { return m+1 ; }
    } ;

  }*/

  template <typename T, typename I=std::ptrdiff_t, typename S=I, int IndexBase=0>
  class crs
  : public crs_structure< T, row_major, I, S, IndexBase >
  {
    private:
      typedef crs_structure< T, row_major, I, S, IndexBase > base_type ;

    public:
      typedef typename base_type::value_type value_type ;
      typedef typename base_type::index_type index_type ;
      typedef typename base_type::size_type  size_type ;

      static size_type const index_base = base_type::index_base ;

    private:
      typedef shared_vector<value_type,S> full_data_type ;
      typedef shared_vector<size_type,S> full_rows_type ;
      typedef shared_vector<index_type,S> full_columns_type ;

    public:
      // Construct CRS from sizes, has empty structure and data vectors
      crs( size_type n_rows, size_type n_columns, size_type nnz )
      : base_type()
      , data_( nnz )
      , rows_( n_rows+1 )
      , columns_( nnz )
      {
        this->reset( n_rows, n_columns, nnz, data_, rows_, columns_ ) ;
      }

    public:
      crs( crs const& that )
      : base_type( that )
      , data_( that.data_ )
      , rows_( that.rows_ )
      , columns_( that.columns_ )
      {}

      template <typename COO>
      crs( COO const& that )
      : base_type()
      , data_( that.num_nz() )
      , rows_( that.num_rows()+1)
      , columns_( that.num_nz())
      {
        assert( is_sorted(that) ) ;
        static_assert( is<CoordinateSparseMatrix,COO>::value, "COO is not a CoordinateSparseMatrix" ) ;

        this->reset( base_type( that.num_rows(), that.num_columns(), that.num_nz(), data_, rows_, columns_) ) ;

        columns_ = that.column_indices() ;
        data_ = that.data() ;

        rows_(0) = 0 ;
        size_type current_row = 0 ;
        for (size_type i=0; i<this->num_nz(); ++i) {
          while (current_row<that.row(i)) {
            current_row++ ;
            rows_(current_row) = i ;
          }
        }
        rows_(this->num_rows()) = this->num_nz() ;
      }

      void resize( size_type n_rows, size_type n_columns, size_type nnz ) {
        data_.resize( nnz ) ;
        rows_.resize( n_rows+1 ) ;
        columns_.resize( nnz ) ;
        this->reset( base_type( n_rows, n_columns, nnz, data_, rows_, columns_ ) ) ;
      }

    private:
      full_data_type    data_ ;
      full_rows_type    rows_ ;
      full_columns_type columns_ ;
  } ;


  template <typename T, typename I, typename S, int IndexBase>
  struct glas_concept< crs<T,I,S,IndexBase> >
  : CompressedSparseMatrix
  {};

} // namespace glas2


#endif
