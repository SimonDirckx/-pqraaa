//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_sparse_container_crs_void_hpp
#define glas2_sparse_container_crs_void_hpp

#include <glas2/sparse/container/crs.hpp>
#include <glas2/sparse/type/crs_structure_void.hpp>

namespace glas2 {

  template <typename I, typename S, int IndexBase>
  class crs< no_value_type, I, S, IndexBase >
  : public crs_structure< no_value_type, row_major, I, S, IndexBase >
  {
    private:
      typedef crs_structure< no_value_type, row_major, I, S, IndexBase > base_type ;

    public:
      typedef typename base_type::value_type value_type ;
      typedef typename base_type::index_type index_type ;
      typedef typename base_type::size_type  size_type ;

      static size_type const index_base = base_type::index_base ;

    private:
      typedef shared_vector<size_type,S> full_rows_type ;
      typedef shared_vector<index_type,S> full_columns_type ;

    public:
      // Construct CRS from sizes, has empty structure and data vectors
      crs( size_type n_rows, size_type n_columns, size_type nnz )
      : base_type()
      , rows_( n_rows+1 )
      , columns_( nnz )
      {
        this->reset( base_type( n_rows, n_columns, nnz, no_value_type_array(), rows_, columns_ ) ) ;
      }

    public:
      crs( crs const& that )
      : base_type( that )
      , rows_( that.rows_ )
      , columns_( that.columns_ )
      {}

      template <typename COO>
      crs( COO const& that )
      : base_type()
      , rows_( that.num_rows()+1)
      , columns_( that.num_nz())
      {
        assert( is_sorted(that) ) ;
        static_assert( is<CoordinateSparseMatrix,COO>::value, "COO is not a CoordinateSparseMatrix" ) ;

        this->reset( base_type( that.num_rows(), that.num_columns(), that.num_nz(), no_value_type_array(), rows_, columns_) ) ;

        columns_ = that.column_indices() ;

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
        rows_.resize( n_rows+1 ) ;
        columns_.resize( nnz ) ;
        this->reset( base_type( n_rows, n_columns, nnz, no_value_type_array(), rows_, columns_ ) ) ;
      }

    private:
      full_rows_type    rows_ ;
      full_columns_type columns_ ;
  } ;


} // namespace glas2


#endif
