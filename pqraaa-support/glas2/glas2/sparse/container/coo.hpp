//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_sparse_container_coo_hpp
#define glas2_sparse_container_coo_hpp

#include <glas2/type/pass_reference.hpp>
#include <glas2/vector/container/shared_vector.hpp>
#include <glas2/sparse/type/coo_structure.hpp>
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
  class coo
  : public coo_structure<T,I,S,IndexBase>
  {
    public:
      typedef coo_structure<T,I,S,IndexBase> base_type ;

    public:
      typedef typename base_type::value_type value_type ;
      typedef typename base_type::index_type index_type ; // E.g. short int
      typedef typename base_type::size_type  size_type ;

    private:
      typedef shared_vector<value_type,S> full_data_type ;
      typedef shared_vector<index_type,S> full_rows_type ;
      typedef shared_vector<index_type,S> full_columns_type ;

    public:
      coo()
      : base_type()
      , full_data_(0)
      , full_rows_(0)
      , full_columns_(0)
      {
        static_cast<base_type&>(*this).reset( 0, 0, 0, full_data_, full_rows_, full_columns_ ) ;
      }

      coo( size_type n, size_type m, size_type reserve_nnz=0 )
      : full_data_(reserve_nnz)
      , full_rows_(reserve_nnz)
      , full_columns_(reserve_nnz)
      {
        static_cast<base_type&>(*this).reset( n, m, 0, full_data_, full_rows_, full_columns_ ) ;
      }

    public: // Copies
      // Shallow copy !!
      coo( coo const& that )
      : base_type(that)
      , full_data_(that.full_data_)
      , full_rows_(that.full_rows_)
      , full_columns_(that.full_columns_)
      {}

    protected:
      void reserve( size_type nnz ) {
        if (nnz>full_data_.size()) {
          full_data_type data_temp( this->num_nz() ) ; data_temp = this->data() ;
          full_data_.resize( nnz ) ; if (this->num_nz()>0) full_data_( range(0,this->num_nz()) ) = data_temp( range(0,this->num_nz()) ) ;
        }
        if (nnz>full_rows_.size()) {
          full_rows_type rows_temp( this->num_nz() ) ; rows_temp = this->row_indices() ;
          full_rows_.resize( nnz ) ; if (this->num_nz()>0) full_rows_( range(0,this->num_nz()) ) = rows_temp( range(0,this->num_nz()) ) ;
        }
        if (nnz>full_columns_.size()) {
          full_columns_type columns_temp( this->num_nz() ) ; columns_temp = this->column_indices() ;
          full_columns_.resize( nnz ) ; if (this->num_nz()>0) full_columns_( range(0,this->num_nz()) ) = columns_temp( range(0,this->num_nz()) ) ;
        }
        static_cast<base_type&>(*this).reset( this->num_rows(), this->num_columns(), this->num_nz(), full_data_, full_rows_, full_columns_ ) ;
      } // reserve()

    public:
      coo& reset( size_type n, size_type m, size_type reserve_nnz=0 ) {
         full_data_.resize( reserve_nnz ) ;
         full_rows_.resize( reserve_nnz ) ;
         full_columns_.resize( reserve_nnz ) ;
         static_cast<base_type&>(*this).reset( n, m, 0, full_data_, full_rows_, full_columns_ ) ;
         return *this ;
      }

      coo& push_back( size_type row, size_type column, value_type value ) {
        if ( this->nnz_ >= full_rows_.size() ) {
          reserve( std::max( this->nnz_ * 2, this->nnz_+this->num_rows()+this->num_columns() ) ) ;
        }
        full_rows_(this->nnz_) = row+IndexBase ;
        full_columns_(this->nnz_) = column+IndexBase ;
        full_data_(this->nnz_) = value ;
        ++this->nnz_ ;
        return *this ;
      }

      template <typename E>
      coo& operator=( E const& that ) {
        return assign( *this, that ) ;
      }

      // we have to override this because operator+ passes the matrix by value
      /*template <typename E>
      coo& operator+=( E const& that ) {
        return plus_assign( current_backend(), *this, that ) ;
      }*/

    private:
      full_data_type            full_data_ ;
      full_rows_type            full_rows_ ;
      full_columns_type         full_columns_ ;
  } ;


  template <typename T, typename I, typename S, int IndexBase>
  struct glas_concept< coo<T,I,S,IndexBase> >
  : ForwardCoordinateSparseMatrix
  {};

  template <typename T, typename I, typename S, int IndexBase>
  struct pass_reference< coo<T,I,S,IndexBase> >
  {
    static const bool pass_by_value = false ;

    typedef coo<T,I,S,IndexBase>& type ;

    type operator()( type that ) const {
      return that ;
    }
  };

} // namespace glas2


#endif
