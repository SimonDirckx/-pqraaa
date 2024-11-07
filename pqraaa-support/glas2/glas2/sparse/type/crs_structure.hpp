//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_sparse_type_crs_structure_hpp
#define glas2_sparse_type_crs_structure_hpp

#include <glas2/sparse/algorithm/assign.hpp>
#include <glas2/sparse/algorithm/is_sorted.hpp>
#include <glas2/sparse/view/sparse_vector_view.hpp>
#include <glas2/sparse/concept/compressed_sparse_matrix.hpp>
#include <glas2/sparse/concept/coordinate_sparse_matrix.hpp>
#include <glas2/matrix/view/matrix_selection.hpp>
#include <glas2/vector/type/contiguous_vector.hpp>
#include <glas2/vector/algorithm/vector_selection.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename T, typename Orientation, typename I=std::ptrdiff_t, typename S=I, int IndexBase=0>
  class crs_structure {
    public:
      typedef T value_type ;
      typedef I index_type ;
      typedef S size_type ; // E.g. short int
      typedef Orientation orientation ;

      static index_type const index_base = IndexBase ;
      typedef contiguous_vector< value_type, size_type > data_type ;
      typedef contiguous_vector< size_type, size_type >  compressed_index_type ;
      typedef contiguous_vector< index_type, size_type > direct_index_type ;

    public:
      crs_structure()
      : num_rows_(0)
      , num_columns_(0)
      , nnz_(0)
      , data_()
      , compressed_index_()
      , direct_index_()
      {}

      crs_structure( compressed_index_type const& compressed_index, direct_index_type const& direct_index, data_type const& data )
      : num_rows_(compressed_index.size()-1)
      , num_columns_(num_rows_)
      , nnz_(direct_index.size())
      , data_(data)
      , compressed_index_( compressed_index )
      , direct_index_( direct_index )
      {
        assert( data.size()==nnz_ ) ;
      }

      crs_structure( size_type n, size_type m, size_type nnz, data_type const& data, compressed_index_type const& compressed_index, direct_index_type const& direct_index )
      : num_rows_(n)
      , num_columns_(m)
      , nnz_(nnz)
      , data_(data)
      , compressed_index_( compressed_index )
      , direct_index_( direct_index )
      {}

      crs_structure( size_type n, size_type m, size_type nnz, value_type* data, size_type* compressed_index, index_type* direct_index )
      : num_rows_(n)
      , num_columns_(m)
      , nnz_(nnz)
      , data_(data)
      , compressed_index_( compressed_index )
      , direct_index_( direct_index )
      {}

      crs_structure( crs_structure const& that )
      : num_rows_( that.num_rows_ )
      , num_columns_( that.num_columns_ )
      , nnz_( that.nnz_ )
      , data_( that.data_ )
      , compressed_index_( that.compressed_index_ )
      , direct_index_( that.direct_index_ )
      {}

      // Matrix
      size_type num_rows() const { return num_rows_ ; }
      size_type num_columns() const { return num_columns_ ; }
      size_type num_nz() const { return nnz_ ; }

      crs_structure& reset( crs_structure const& that ) {
        num_rows_ = that.num_rows_ ;
        num_columns_ = that.num_columns_ ;
        nnz_ = that.nnz_ ;
        data_.reset( that.data_.ptr(), that.data_.size() ) ;
        compressed_index_.reset( that.compressed_index_.ptr(), that.compressed_index_.size() ) ;
        direct_index_.reset( that.direct_index_.ptr(), that.direct_index_.size() ) ;
        return *this ;
      }

      crs_structure& reset( size_type n, size_type m, size_type nnz, data_type const& data, compressed_index_type const& compressed_index, direct_index_type const& direct_index ) {
        num_rows_ = n ;
        num_columns_ = m ;
        nnz_ = nnz ;
        data_.reset( data.ptr(), data.size() ) ;
        compressed_index_.reset( compressed_index.ptr(), compressed_index.size() ) ;
        direct_index_.reset( direct_index.ptr(), direct_index.size() ) ;
        return *this ;
      }

      crs_structure& operator=( crs_structure const& that ) {
        assert( that.num_rows()==num_rows() ) ;
        assert( that.num_columns()==num_columns() ) ;
        assert( that.num_nz()==num_nz() ) ;
        data_ = that.data_ ;
        compressed_index_ = that.compressed_index_ ;
        direct_index_ = that.direct_index_ ;
        return *this ;
      }

      template <typename E>
      crs_structure operator=( E const& that ) {
        return assign( *this, that ) ;
      }

      data_type data() { return data_ ; }
      data_type data() const { return data_ ; }
      compressed_index_type compressed_index() const { return compressed_index_ ; }
      direct_index_type direct_index() const { return direct_index_ ; }

    public:
      template <typename I1, typename I2>
      typename std::enable_if< !(std::is_integral<I1>::value && std::is_integral<I2>::value)
                             , typename matrix_selection< crs_structure, I1, I2 >::result_type
                             >::type operator()( I1 const& s1, I2 const& s2 ) const {
        return matrix_selection< crs_structure, I1, I2 >::apply( *this, s1, s2 ) ;
      }

    private:
      size_type             num_rows_ ;
      size_type             num_columns_ ;
      size_type             nnz_ ;
      data_type             data_ ;
      compressed_index_type compressed_index_ ;
      direct_index_type    direct_index_ ;
  } ;


  template <typename T, typename O, typename I, typename S, int IndexBase>
  struct glas_concept< crs_structure<T,O,I,S,IndexBase> >
  : CompressedSparseMatrix
  {};


  template <typename T, typename I, typename S, int IndexBase, typename Row>
  struct matrix_selection< crs_structure<T,row_major,I,S,IndexBase>, Row, all,
        typename std::enable_if< std::is_integral<Row>::value >::type > {
    typedef typename vector_selection< typename crs_structure<T,row_major,I,S,IndexBase>::data_type, range >::result_type         data_type ;
    typedef typename vector_selection< typename crs_structure<T,row_major,I,S,IndexBase>::direct_index_type, range >::result_type index_type ;
    typedef sparse_vector_view< data_type, index_type, IndexBase > result_type ;

    static result_type apply( crs_structure<T,row_major,I,S,IndexBase> m, Row const& r, all ) {
      assert( r<m.num_rows() && r>=0 ) ;
      range r_row( m.compressed_index()(r), m.compressed_index()(r+1) ) ;
      return result_type( m.data()( r_row ), m.direct_index()( r_row ) , m.num_columns() ) ;
    }
  } ;

} // namespace glas2


#endif
