//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_sparse_type_coo_structure_hpp
#define glas2_sparse_type_coo_structure_hpp

#include <glas2/vector/type/contiguous_vector.hpp>
#include <glas2/scalar/container/shared_scalar.hpp>
#include <glas2/sparse/algorithm/assign.hpp>
#include <glas2/sparse/algorithm/ops_assign.hpp>
#include <glas2/sparse/algorithm/sort_and_compress.hpp>
#include <glas2/sparse/concept/forward_coordinate_sparse_matrix.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename T, typename I=std::ptrdiff_t, typename S=I, int IndexBase=0>
  class coo_structure {
    public:
      typedef T value_type ;
      typedef I index_type ;
      typedef S size_type ; // E.g. short int

      static index_type const index_base = IndexBase ;

      typedef contiguous_vector< value_type, size_type > contiguous_data_type ;
      typedef contiguous_vector< index_type, size_type >  contiguous_rows_type ;
      typedef contiguous_vector< index_type, size_type > contiguous_columns_type ;

    public:
      coo_structure()
      : num_rows_(0)
      , num_columns_(0)
      , nnz_(0)
      , data_()
      , rows_()
      , columns_()
      {}

      coo_structure( size_type n, contiguous_rows_type const& rows, contiguous_columns_type const& columns, contiguous_data_type const& data )
      : num_rows_(n)
      , num_columns_(n)
      , nnz_(data.size())
      , data_(data)
      , rows_( rows )
      , columns_( columns )
      {
        assert( nnz_==rows.size() ) ;
        assert( nnz_==columns.size() ) ;
      }

      coo_structure( size_type n, size_type m, size_type nnz, contiguous_data_type const& data, contiguous_rows_type const& rows, contiguous_columns_type const& columns )
      : num_rows_(n)
      , num_columns_(m)
      , nnz_(nnz)
      , data_(data)
      , rows_( rows )
      , columns_( columns )
      {
        assert( nnz<=rows.size() ) ;
        assert( nnz<=columns.size() ) ;
        assert( nnz<=data.size() ) ;
      }

      coo_structure( size_type n, size_type m, size_type nnz, value_type* data, index_type* rows, index_type* columns )
      : num_rows_(n)
      , num_columns_(m)
      , nnz_(nnz)
      , data_(data,nnz)
      , rows_( rows,nnz )
      , columns_( columns,nnz )
      {}

    public: // Copies
      // Shallow copy !!
      coo_structure( coo_structure const& that )
      : num_rows_(that.num_rows_)
      , num_columns_(that.num_columns_)
      , nnz_(that.nnz_)
      , data_(that.data_)
      , rows_(that.rows_)
      , columns_(that.columns_)
      {}

    public:
      // Matrix
      size_type num_rows() const { return num_rows_ ; }
      size_type num_columns() const { return num_columns_ ; }
      size_type num_nz() const { return nnz_ ; }

      coo_structure& operator=( coo_structure const& that ) {
        assert( that.num_rows()==num_rows() ) ;
        assert( that.num_columns()==num_columns() ) ;
        assert( that.num_nz()==num_nz() ) ;
        data_ = that.data_ ;
        rows_ = that.rows_ ;
        columns_ = that.columns_ ;
        return *this ;
      }

      coo_structure& reset( size_type n, size_type m, size_type nnz, contiguous_data_type data, contiguous_rows_type rows, contiguous_columns_type columns ) {
        num_rows_ = n ;
        num_columns_ = m ;
        nnz_ = nnz ;
        data_.reset( data.ptr(), data.size() ) ;
        rows_.reset( rows.ptr(), rows.size() ) ;
        columns_.reset( rows.ptr(), rows.size() ) ;
        columns_.reset( columns.ptr(), columns.size() ) ;
        assert( nnz<=rows.size() ) ;
        assert( nnz<=columns.size() ) ;
        assert( nnz<=data.size() ) ;
        return *this ;
      }

      typedef typename vector_selection<contiguous_data_type,range>::result_type data_type ;
      data_type data() const { return data_(range(0,nnz_) ) ; }
      data_type data() { return data_(range(0,nnz_) ) ; }

      typedef typename vector_selection<contiguous_rows_type,range>::result_type rows_type ;
      rows_type row_indices() const { return rows_(range(0,nnz_) ) ; }

      typedef typename vector_selection<contiguous_columns_type,range>::result_type columns_type ;
      columns_type column_indices() const { return columns_(range(0,nnz_) ) ; }

      size_type row( size_type i ) const {
        assert( i>=0 && i<nnz_ ) ;
        return rows_(i) - index_base ;
      }

      size_type column( size_type i ) const {
        assert( i>=0 && i<nnz_ ) ;
        return columns_(i) - index_base ;
      }

      coo_structure& sort_and_compress() {
        glas2::sort_and_compress( *this ) ;

        // Merge doubles
        size_type offset = 1 ;
        for (size_type i=offset ; i<nnz_ ; ++i) {
          if (row(i)==row(i-offset) && column(i)==column(i-offset)) {
            data()(i-offset) += data()(i) ;
            ++offset ;
          } else {
            data_(i-offset+1) = data_(i) ;
            rows_(i-offset+1) = rows_(i) ;
            columns_(i-offset+1) = columns_(i) ;
          }
        }
        nnz_ -= offset - 1 ;
        return *this ;
      }

    public:
      class iterator {
        public:
          typedef typename coo_structure::value_type value_type ;
          typedef typename coo_structure::size_type  size_type ;

          iterator( size_type index, contiguous_rows_type const& rows, contiguous_columns_type const& columns, contiguous_data_type const& data )
          : index_( index )
          , rows_( rows )
          , columns_( columns )
          , data_( data )
          {}

        public:
          void operator++() { ++index_ ; }
          size_type row() const { return rows_(index_)-index_base ; }
          size_type column() const { return columns_(index_)-index_base ; }
          value_type value() const { return data_(index_) ; }

          bool operator==( iterator const& that ) const { return index_ == that.index_ ; }
          bool operator!=( iterator const& that ) const { return index_ != that.index_ ; }

        private:
          size_type           index_ ;
          contiguous_rows_type const&    rows_ ;
          contiguous_columns_type const& columns_ ;
          contiguous_data_type const&    data_ ;
      } ;

      iterator begin() const { return iterator( 0, rows_, columns_, data_ ) ; }
      iterator end() const { return iterator( nnz_, rows_, columns_, data_ ) ; }

    public:
      contiguous_data_type const& storage_data() const { return data_ ; }
      contiguous_rows_type const& storage_rows() const { return rows_ ; }
      contiguous_columns_type const& storage_columns() const { return columns_ ; }

    private:
      size_type                 num_rows_ ;
      size_type                 num_columns_ ;

    protected:
      shared_scalar<size_type>  nnz_ ;
      contiguous_data_type      data_ ;
      contiguous_rows_type      rows_ ;
      contiguous_columns_type   columns_ ;
  } ;


  template <typename T, typename I, typename S, int IndexBase>
  struct glas_concept< coo_structure<T,I,S,IndexBase> >
  : ForwardCoordinateSparseMatrix
  {};

} // namespace glas2


#endif
