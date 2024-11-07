//  (C) Copyright Karl Meerbergen 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_sparse_container_coo_void_hpp
#define glas2_sparse_container_coo_void_hpp

#include <glas2/sparse/type/no_value_type_array.hpp>
#include <glas2/sparse/container/coo.hpp>

namespace glas2 {

  template <typename I, typename S, int IndexBase>
  class coo< no_value_type, I, S, IndexBase > {
    public:
      typedef no_value_type value_type ;
      typedef I             index_type ;
      typedef S             size_type ; // E.g. short int

      static size_type const index_base = IndexBase ;

    private:
      typedef shared_vector<index_type,S> full_rows_type ;
      typedef shared_vector<index_type,S> full_columns_type ;

    public:
      coo()
      : num_rows_(0)
      , num_columns_(0)
      , nnz_(0)
      , rows_(0)
      , columns_(0)
      {}

      coo( size_type n, size_type m, size_type reserve_nnz=0 )
      : num_rows_(n)
      , num_columns_(m)
      , nnz_(0)
      , rows_(reserve_nnz)
      , columns_(reserve_nnz)
      {}

      // Deep copy !!
      coo( coo const& that )
      : num_rows_(that.num_rows_)
      , num_columns_(that.num_columns_)
      , nnz_(that.nnz_)
      , rows_(that.rows_)
      , columns_(that.columns_)
      {}

    protected:
      void reserve( size_type nnz ) {
        if (nnz>rows_.size()) {
          full_rows_type rows_temp( nnz_ ) ; rows_temp = row_indices() ;
          rows_.resize( nnz ) ; rows_( range(0,nnz_) ) = rows_temp( range(0,nnz_) ) ;
        }
        if (nnz>columns_.size()) {
          full_columns_type columns_temp( nnz_ ) ; columns_temp = column_indices() ;
          columns_.resize( nnz ) ; columns_( range(0,nnz_) ) = columns_temp( range(0,nnz_) ) ;
        }
      } // reserve()

    public:
      coo& reset( size_type n, size_type m, size_type reserve_nnz=0 ) {
         num_rows_ = n ;
         num_columns_ = m ;
         rows_.resize( reserve_nnz ) ;
         columns_.resize( reserve_nnz ) ;
         nnz_ = 0 ;
      }

      // Matrix
      size_type num_rows() const { return num_rows_ ; }
      size_type num_columns() const { return num_columns_ ; }
      size_type num_nz() const { return nnz_ ; }

      coo& operator=( coo const& that ) {
        assert( that.num_rows()==num_rows() ) ;
        assert( that.num_columns()==num_columns() ) ;
        assert( that.num_nz()==num_nz() ) ;
        rows_ = that.rows_ ;
        columns_ = that.columns_ ;
        return *this ;
      }

      template <typename E>
      coo& operator=( E const& that ) {
        return assign( *this, that ) ;
      }

      template <typename E>
      coo& operator+=( E const& that ) {
        return plus_assign( current_backend(), *this, that ) ;
      }

      coo& push_back( size_type row, size_type column, value_type v=value_type() ) {
        if ( nnz_ >= rows_.size() ) {
          reserve( std::max( nnz_ * 2, nnz_+num_rows_+num_columns_ ) ) ;
        }
        rows_(nnz_) = row+IndexBase ;
        columns_(nnz_) = column+IndexBase ;
        ++nnz_ ;
        return *this ;
      }

      typedef no_value_type_array data_type ;
      data_type data() const { return no_value_type_array() ; }

      typedef typename vector_selection<typename full_rows_type::base_type,range>::result_type rows_type ;
      rows_type row_indices() const { return rows_(range(0,nnz_) ) ; }

      typedef typename vector_selection<typename full_columns_type::base_type,range>::result_type columns_type ;
      columns_type column_indices() const { return columns_(range(0,nnz_) ) ; }

      size_type row( size_type i ) const {
        assert( i>=0 && i<nnz_ ) ;
        return rows_(i) - index_base ;
      }

      size_type column( size_type i ) const {
        assert( i>=0 && i<nnz_ ) ;
        return columns_(i) - index_base ;
      }

      coo& sort_and_compress() {
        // Bubble sort
        bool done = false ;

        while (!done) {
          size_type offset = 1 ;
          done = true ;
          for (size_type i=offset ; i<nnz_ ; ++i) {
            if (offset>1) {
              rows_(i-offset+1) = rows_(i) ;
              columns_(i-offset+1) = columns_(i) ;
            }
            if (rows_(i-offset) > rows_(i-offset+1) ||
                 ( rows_(i-offset) == rows_(i-offset+1) && columns_(i-offset) > columns_(i-offset+1) )
               ) {
              done = false ;
              size_type ii = rows_(i-offset+1) ; rows_(i-offset+1) = rows_(i-offset) ; rows_(i-offset) = ii ;
              ii = columns_(i-offset+1) ; columns_(i-offset+1) = columns_(i-offset) ; columns_(i-offset) = ii ;
            } else if (rows_(i-offset) == rows_(i-offset+1) && columns_(i-offset) == columns_(i-offset+1) ) {
              offset ++ ;
            }
          }
          nnz_ -= (offset-1) ;
        }
        return *this ;
      }

    private:
      size_type         num_rows_ ;
      size_type         num_columns_ ;
      size_type         nnz_ ;
      full_rows_type    rows_ ;
      full_columns_type columns_ ;
  } ;

} // namespace glas2


#endif
