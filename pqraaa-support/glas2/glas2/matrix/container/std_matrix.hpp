//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#include <glas2/matrix/algorithm/assign.hpp>
#include <glas2/matrix/iterator/contiguous_iterator.hpp>
#include <glas2/matrix/view/matrix_selection.hpp>
#include <glas2/matrix/concept/contiguous_dense_matrix.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/concept/is.hpp>
#include <algorithm>
#include <type_traits>
#include <cassert>

namespace glas2 {

  namespace detail_std_matrix {
    template <typename S>
    void select_contiguous_iterator( row_major, S& size_1, S& size_2, S num_rows_, S num_columns_ ) {
      size_1 = num_rows_ ;
      size_2 = num_columns_ ;
    }

    template <typename S>
    void select_contiguous_iterator( column_major, S& size_1, S& size_2, S num_rows_, S num_columns_ ) {
      size_1 = num_columns_ ;
      size_2 = num_rows_ ;
    }

    template <typename S>
    S address( row_major, S i, S j, S nr, S nc ) { return i*nc+j ; }

    template <typename S>
    S address( column_major, S i, S j, S nr, S nc ) { return j*nr+i ; }

  } // namespace detail_std_matrix

  template <typename V, typename O=column_major>
  class std_matrix {
    public:
      typedef V                                vector_type ;
      typedef O                                orientation ;
      typedef general_matrix                   structure ;
      typedef typename vector_type::value_type value_type ;
      typedef typename vector_type::size_type  size_type ;

      typedef ContiguousDenseMatrix            concept ;

    public:
        std_matrix( vector_type& v, size_type num_rows, size_type num_columns )
        : v_(v)
        , num_rows_(num_rows)
        , num_columns_(num_columns)
        {
          assert( v_.size() == num_columns*num_rows ) ;
        }

        // Copy reference !!
        std_matrix( std_matrix const& that )
        : v_(that.v_)
        , num_rows_(that.num_rows_)
        , num_columns_(that.num_columns_)
        {}


    public:
      size_type num_rows() const { return num_rows_ ; }
      size_type num_columns() const { return num_columns_ ; }

      value_type const& operator() ( size_type i, size_type j ) const { return v_[detail_std_matrix::address(orientation(),i,j,num_rows_,num_columns_)] ; }
      value_type& operator() ( size_type i, size_type j ) { return v_[detail_std_matrix::address(orientation(),i,j,num_rows_,num_columns_)] ; }

      size_type stride() const { return 1 ; }
      value_type* storage_ptr() { return &v_.begin() ; }
      value_type const* storage_ptr() const { return &v_.begin() ; }

      std_matrix<V>& operator=( std_matrix const& that ) {
        assert( that.v_.size()()==v_.size() ) ;
        std::copy( that.storage_ptr(), that.storage_ptr()+that.v_.size(), storage_ptr() ) ;
        return *this ;
      }

      template <typename E>
      std_matrix<V> operator=( E const& that ) {
        return assign( *this, that ) ;
      }

      template <typename S>
      typename std::enable_if< is<DenseMatrix,S>::value, matrix_selection< std_matrix<V>, S1, S2 > >::type operator()( S1 const& s1, S2 const& s2 ) {
        return matrix_selection< std_matrix<V>, S1, S2 >( *this, s1, s2 ) ;
      }

    public:
      typedef contiguous_iterator<size_type,orientation> iterator ;
      iterator iterate() const {
        size_type size_1, size_2 ;
        detail_std_matrix::select_contiguous_iterator( orientation(), size_1, size_2, num_rows_, num_columns_ ) ;
        return contiguous_iterator<size_type,orientation>( size_1, size_2 ) ;
      }

      value_type operator[](size_type i) const { return v_[i]; }
      value_type& operator[](size_type i) { return v_[i]; }

    private:
      vector_type& v_ ;
      size_type    num_rows_ ;
      size_type    num_columns_ ;
  } ;


  template <typename V, typename O>
  struct glas_concept< std_matrix<V,O> >
  : ContiguousDenseMatrix
  {};

} // namespace glas2
