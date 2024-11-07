//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_sparse_type_crs_structure_void_hpp
#define glas2_sparse_type_crs_structure_void_hpp

#include <glas2/sparse/type/no_value_type_array.hpp>
#include <glas2/sparse/type/crs_structure.hpp>
#include <glas2/sparse/view/sparse_vector_view_void.hpp>

namespace glas2 {

  template <typename Orientation, typename I, typename S, int IndexBase>
  class crs_structure< no_value_type, Orientation,I,S,IndexBase> {
    public:
      typedef no_value_type value_type ;
      typedef I             index_type ;
      typedef S             size_type ;
      typedef Orientation   orientation ;

      static index_type const index_base = IndexBase ;
      typedef no_value_type_array                        data_type ;
      typedef contiguous_vector< size_type, size_type >  compressed_index_type ;
      typedef contiguous_vector< index_type, size_type > direct_index_type ;

    public:
      crs_structure()
      : num_rows_(0)
      , num_columns_(0)
      , nnz_(0)
      , compressed_index_()
      , direct_index_()
      {}

      crs_structure( size_type n, size_type m, size_type nnz, data_type const& data, compressed_index_type const& compressed_index, direct_index_type const& direct_index )
      : num_rows_(n)
      , num_columns_(m)
      , nnz_(nnz)
      , compressed_index_( compressed_index )
      , direct_index_( direct_index )
      {}

      crs_structure( crs_structure const& that )
      : num_rows_( that.num_rows_ )
      , num_columns_( that.num_columns_ )
      , nnz_( that.nnz_ )
      , compressed_index_( that.compressed_index_ )
      , direct_index_( that.direct_index_ )
      {}

    public:
      // Matrix
      size_type num_rows() const { return num_rows_ ; }
      size_type num_columns() const { return num_columns_ ; }
      size_type num_nz() const { return nnz_ ; }

      crs_structure& reset( crs_structure const& that ) {
        num_rows_ = that.num_rows_ ;
        num_columns_ = that.num_columns_ ;
        nnz_ = that.nnz_ ;
        compressed_index_.reset( that.compressed_index_.ptr(), that.compressed_index_.size() ) ;
        direct_index_.reset( that.direct_index_.ptr(), that.direct_index_.size() ) ;
        return *this ;
      }

      crs_structure& operator=( crs_structure const& that ) {
        assert( that.num_rows()==num_rows() ) ;
        assert( that.num_columns()==num_columns() ) ;
        assert( that.num_nz()==num_nz() ) ;
        compressed_index_ = that.compressed_index_ ;
        direct_index_ = that.direct_index_ ;
        return *this ;
      }

      template <typename E>
      crs_structure operator=( E const& that ) {
        return assign( *this, that ) ;
      }

      data_type data() const { return data_type() ; }
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
      compressed_index_type compressed_index_ ;
      direct_index_type    direct_index_ ;
  } ;


  template <typename I, typename S, int IndexBase, typename Row>
  struct matrix_selection< crs_structure<no_value_type,row_major,I,S,IndexBase>, Row, all,
        typename std::enable_if< std::is_integral<Row>::value >::type > {
    typedef typename vector_selection< typename crs_structure<no_value_type,row_major,I,S,IndexBase>::direct_index_type, range >::result_type index_type ;
    typedef sparse_vector_view< no_value_type_array, index_type, IndexBase > result_type ;

    static result_type apply( crs_structure<no_value_type,row_major,I,S,IndexBase> m, Row const& r, all ) {
      assert( r<m.num_rows() && r>=0 ) ;
      range r_row( m.compressed_index()(r), m.compressed_index()(r+1) ) ;
      return result_type( m.direct_index()( r_row ) , m.num_columns() ) ;
    }
  } ;



} // namespace glas2


#endif
